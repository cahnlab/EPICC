"""
Post-run validation tests for S. pombe pipeline outputs.

These tests validate that a completed pipeline run produced the expected
files with correct formats. They require an existing completed run in
results/ and are skipped automatically if no results are found.

Run with: pytest tests/integration/test_pombe_postrun.py -v -m slow
"""

import gzip
import json
import struct
import sys

import pytest
import yaml
from pathlib import Path

from tests.integration.conftest import DATA_DIR, REPO_ROOT, load_output_dir

sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))
from sample_sheet import (  # noqa: E402
    add_compat_columns,
    build_analysis_name,
    get_analysis_samples,
    get_replicate_sample_ids,
    is_peak_call_target,
    read_sample_sheet,
)

# ---------------------------------------------------------------------------
# Constants: repo root, results directory, run identity
# ---------------------------------------------------------------------------
OPTIONS_FILE = "test_options_pombe.yaml"
RESULTS = REPO_ROOT / load_output_dir(OPTIONS_FILE)

with open(DATA_DIR / OPTIONS_FILE) as _fh:
    _OPTIONS = yaml.safe_load(_fh)
ANALYSIS_NAME = _OPTIONS["analysis_name"]

# ---------------------------------------------------------------------------
# Expected names, derived from the sample sheet with the pipeline's own helpers
# ---------------------------------------------------------------------------
# Naming splits at the alignment boundary (issue #39): download/trim outputs
# use the bare Sample_ID, but every per-replicate file from alignment onward is
# named by mapped_name = '{Sample_ID}__{Genome}' (final__*.bam, FC__final__*.bw,
# peaks_*__final__*, map_<env>__*.done, RNA/sRNA tracks, ShortStack dirs).
# Analysis-level (merged) files use build_analysis_name, which already ends in
# the genome. Building both from the sheet keeps this test in step with the
# rules rather than with a hand-copied list (issue #81).
SAMPLES = add_compat_columns(
    read_sample_sheet(REPO_ROOT / _OPTIONS["sample_file"])
)
_GENOMES = sorted(SAMPLES["ref_genome"].unique())
assert len(_GENOMES) == 1, f"pombe postrun expects one genome, got {_GENOMES}"
GENOME = _GENOMES[0]


def _env_rows(env):
    return SAMPLES[SAMPLES["env"] == env]


def _mapped_names(rows):
    return rows["mapped_name"].tolist()


def _merged_analysis_rows(env):
    """Analysis names with >=2 replicates, i.e. those that get merged outputs.

    A single-replicate group has no merged BAM: its downstream files fall back
    to the replicate's own name (see combined_analysis.smk).
    """
    analysis = add_compat_columns(get_analysis_samples(_env_rows(env)))
    rows = {build_analysis_name(row): row for _, row in analysis.iterrows()}
    return {n: r for n, r in rows.items()
            if len(get_replicate_sample_ids(n, SAMPLES)) >= 2}


def _merged_analysis_names(env):
    return list(_merged_analysis_rows(env))


_CHIP_ROWS = _env_rows("ChIP")
_CHIP_PEAK_ROWS = _CHIP_ROWS[_CHIP_ROWS.apply(is_peak_call_target, axis=1)]

# ChIP per-replicate BAMs: every ChIP row, controls included
CHIP_PER_REP_BAMS = _mapped_names(_CHIP_ROWS)

# ChIP per-replicate rows that call peaks (have a Control) -> FC bigwigs
CHIP_NON_INPUT_PER_REP = _mapped_names(_CHIP_PEAK_ROWS)

# Per-replicate peak files: PE -> peaks_pe__*.broadPeak, SE -> peaks_se__*.narrowPeak
CHIP_PE_PEAK_SAMPLES = _mapped_names(
    _CHIP_PEAK_ROWS[_CHIP_PEAK_ROWS["Read_layout"] == "PE"]
)
CHIP_SE_PEAK_SAMPLES = _mapped_names(
    _CHIP_PEAK_ROWS[_CHIP_PEAK_ROWS["Read_layout"] == "SE"]
)
PEAK_TYPE = dict(zip(_CHIP_PEAK_ROWS["mapped_name"], _CHIP_PEAK_ROWS["Peak_type"]))

# ChIP analysis groups with >=2 replicates (merged BAM / FC bigwig / IDR)
CHIP_MERGED_ANALYSIS_BAMS = _merged_analysis_names("ChIP")
MERGED_PEAK_TYPE = {n: r["Peak_type"] for n, r in _merged_analysis_rows("ChIP").items()}

# RNA
RNA_PER_REP = _mapped_names(_env_rows("RNA"))
RNA_MERGED = _merged_analysis_names("RNA")

# sRNA
SRNA_PER_REP = _mapped_names(_env_rows("sRNA"))
SRNA_MERGED = _merged_analysis_names("sRNA")

# sRNA size classes
SRNA_SIZES = ["21nt", "22nt", "23nt", "24nt"]
SRNA_STRANDS = ["plus", "minus"]

# Checkpoint files
CHIP_CHECKPOINTS = [
    f"ChIP_analysis__{ANALYSIS_NAME}__{GENOME}.done",
] + [f"idr__{name}.done" for name in CHIP_MERGED_ANALYSIS_BAMS] \
  + [f"map_ChIP__{s}.done" for s in CHIP_PER_REP_BAMS]

RNA_CHECKPOINTS = [
    f"calling_DEGs__{ANALYSIS_NAME}__{GENOME}.done",
    f"RNA_analysis__{ANALYSIS_NAME}__{GENOME}.done",
] + [f"map_RNA__{s}.done" for s in RNA_PER_REP]

SRNA_CHECKPOINTS = [
    f"sRNA_analysis__{ANALYSIS_NAME}__{GENOME}.done",
    f"calling_differential_sRNA_clusters__{ANALYSIS_NAME}__{GENOME}__on_all_genes.done",
    f"calling_differential_sRNA_clusters__{ANALYSIS_NAME}__{GENOME}__on_new_clusters.done",
] + [f"map_sRNA__{s}.done" for s in SRNA_PER_REP]

# '{analysis_name}__{genome}': the run-level tag on DEG, cluster and combined outputs
RUN = f"{ANALYSIS_NAME}__{GENOME}"

COMBINED_CHECKPOINTS = [
    f"combined_analysis__{ANALYSIS_NAME}__{GENOME}.done",
    f"final_analysis__{ANALYSIS_NAME}.done",
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _file_exists_nonempty(path):
    """Check that a file exists and is non-empty."""
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def _read_bigwig_magic(path):
    """Read the first 4 bytes and return the magic number."""
    with open(path, "rb") as f:
        return struct.unpack("<I", f.read(4))[0]


# ---------------------------------------------------------------------------
# Guard fixture: skip all tests if no completed run exists
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def results_exist():
    """Skip all tests if results/ has no completed pipeline output."""
    chip_dir = RESULTS / "ChIP" / "mapped"
    rna_dir = RESULTS / "RNA" / "mapped"
    if not chip_dir.exists() or not any(chip_dir.glob("final__*.bam")):
        pytest.skip("No completed pipeline run found in results/")
    if not rna_dir.exists() or not any(rna_dir.glob("final__*.bam")):
        pytest.skip("No completed pipeline run found in results/")


# ===========================================================================
# TestChIPOutputs
# ===========================================================================

@pytest.mark.slow
@pytest.mark.integration
class TestChIPOutputs:
    """Validate ChIP-seq pipeline outputs."""

    def test_per_rep_bams_exist(self, results_exist):
        """Every ChIP per-replicate BAM and BAI exists and is non-empty."""
        for sample in CHIP_PER_REP_BAMS:
            bam = RESULTS / "ChIP" / "mapped" / f"final__{sample}.bam"
            bai = RESULTS / "ChIP" / "mapped" / f"final__{sample}.bam.bai"
            assert _file_exists_nonempty(bam), f"Missing or empty: {bam}"
            assert _file_exists_nonempty(bai), f"Missing or empty: {bai}"

    def test_merged_bams_exist(self, results_exist):
        """Merged BAMs are temp() — verify via downstream FC bigwigs instead."""
        for name in CHIP_MERGED_ANALYSIS_BAMS:
            bw = RESULTS / "ChIP" / "tracks" / f"FC__merged__{name}.bw"
            assert _file_exists_nonempty(bw), (
                f"Merged BAM for {name} was never created (downstream bigwig missing): {bw}"
            )

    def test_fc_bigwigs_per_rep(self, results_exist):
        """FC bigwigs exist for every per-rep sample that has a Control."""
        for sample in CHIP_NON_INPUT_PER_REP:
            bw = RESULTS / "ChIP" / "tracks" / f"FC__final__{sample}.bw"
            assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_fc_bigwigs_merged(self, results_exist):
        """Merged FC bigwigs exist for every multi-replicate analysis group."""
        for name in CHIP_MERGED_ANALYSIS_BAMS:
            bw = RESULTS / "ChIP" / "tracks" / f"FC__merged__{name}.bw"
            assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_broad_peaks(self, results_exist):
        """Per-replicate peak files for PE samples (peaks_pe__) exist."""
        for sample in CHIP_PE_PEAK_SAMPLES:
            peak = RESULTS / "ChIP" / "peaks" / f"peaks_pe__final__{sample}_peaks.{PEAK_TYPE[sample]}Peak"
            assert _file_exists_nonempty(peak), f"Missing or empty: {peak}"

    def test_narrow_peaks(self, results_exist):
        """Per-replicate peak files for SE samples (peaks_se__) exist."""
        for sample in CHIP_SE_PEAK_SAMPLES:
            peak = RESULTS / "ChIP" / "peaks" / f"peaks_se__final__{sample}_peaks.{PEAK_TYPE[sample]}Peak"
            assert _file_exists_nonempty(peak), f"Missing or empty: {peak}"

    def test_idr_files(self, results_exist):
        """IDR output files exist (may be empty when IDR fails due to few peaks)."""
        for name in CHIP_MERGED_ANALYSIS_BAMS:
            peaktype = MERGED_PEAK_TYPE[name]
            idr = RESULTS / "ChIP" / "peaks" / f"idr_peaks__{name}.{peaktype}Peak"
            assert idr.exists(), f"Missing: {idr}"

    def test_chip_checkpoints(self, results_exist):
        """All ChIP checkpoints present."""
        for chkpt in CHIP_CHECKPOINTS:
            path = RESULTS / "ChIP" / "chkpts" / chkpt
            assert path.exists(), f"Missing checkpoint: {path}"


# ===========================================================================
# TestRNAOutputs
# ===========================================================================

@pytest.mark.slow
@pytest.mark.integration
class TestRNAOutputs:
    """Validate RNA-seq pipeline outputs."""

    def test_per_rep_bams(self, results_exist):
        """RNA per-replicate BAMs + BAIs exist."""
        for sample in RNA_PER_REP:
            bam = RESULTS / "RNA" / "mapped" / f"final__{sample}.bam"
            bai = RESULTS / "RNA" / "mapped" / f"final__{sample}.bam.bai"
            assert _file_exists_nonempty(bam), f"Missing or empty: {bam}"
            assert _file_exists_nonempty(bai), f"Missing or empty: {bai}"

    def test_merged_bams(self, results_exist):
        """Merged BAMs are temp() — verify via downstream strand bigwigs."""
        for name in RNA_MERGED:
            bw = RESULTS / "RNA" / "tracks" / f"{name}__plus.bw"
            assert _file_exists_nonempty(bw), (
                f"Merged BAM for {name} was never created (downstream bigwig missing): {bw}"
            )

    def test_stranded_bigwigs(self, results_exist):
        """Strand-specific bigwigs for per-rep and merged samples."""
        for sample in RNA_PER_REP:
            for strand in ["plus", "minus"]:
                bw = RESULTS / "RNA" / "tracks" / f"{sample}__{strand}.bw"
                assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"
        for name in RNA_MERGED:
            for strand in ["plus", "minus"]:
                bw = RESULTS / "RNA" / "tracks" / f"{name}__{strand}.bw"
                assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_count_matrix(self, results_exist):
        """RPKM gene expression matrix has correct column structure.

        Note: the raw count matrix (counts__*.txt) is temp() and removed after
        DEG analysis, so we validate the persistent RPKM output instead.
        """
        rpkm = RESULTS / "RNA" / "DEG" / f"genes_rpkm__{RUN}.txt"
        assert _file_exists_nonempty(rpkm), f"Missing or empty: {rpkm}"
        with open(rpkm) as f:
            header = f.readline().strip().split("\t")
        assert header[0] == "GID", f"First column should be GID, got: {header[0]}"
        assert "RPKM" in header, f"RPKM column missing from header: {header}"

    def test_deg_file(self, results_exist):
        """DEG file is non-empty and has expected header columns.

        Note: with a single factor (genotype), sample labels get a trailing
        underscore from the empty 2nd factor, producing WT___vs_dcr1__.
        """
        deg = RESULTS / "RNA" / "DEG" / f"DEG_{RUN}__WT___vs_dcr1__.txt"
        assert _file_exists_nonempty(deg), f"Missing or empty: {deg}"
        with open(deg) as f:
            header = f.readline().strip()
        assert "logFC" in header, f"DEG header missing logFC: {header}"
        assert "DEG" in header, f"DEG header missing DEG column: {header}"

    def test_rna_checkpoints(self, results_exist):
        """All RNA checkpoints present."""
        for chkpt in RNA_CHECKPOINTS:
            path = RESULTS / "RNA" / "chkpts" / chkpt
            assert path.exists(), f"Missing checkpoint: {path}"


# ===========================================================================
# TestSRNAOutputs
# ===========================================================================

@pytest.mark.slow
@pytest.mark.integration
class TestSRNAOutputs:
    """Validate small RNA pipeline outputs."""

    def test_shortstack_results(self, results_exist):
        """ShortStack Results.txt for every sRNA replicate."""
        for sample in SRNA_PER_REP:
            results_txt = RESULTS / "sRNA" / "mapped" / sample / "Results.txt"
            assert _file_exists_nonempty(results_txt), f"Missing or empty: {results_txt}"

    def test_size_class_bigwigs_per_rep(self, results_exist):
        """Size-class bigwigs (4 sizes x 2 strands x every replicate)."""
        for sample in SRNA_PER_REP:
            for size in SRNA_SIZES:
                for strand in SRNA_STRANDS:
                    bw = RESULTS / "sRNA" / "tracks" / f"{sample}__{size}__{strand}.bw"
                    assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_merged_wt_bigwigs(self, results_exist):
        """Merged bigwigs for multi-replicate groups (4 sizes x 2 strands)."""
        for name in SRNA_MERGED:
            for size in SRNA_SIZES:
                for strand in SRNA_STRANDS:
                    bw = RESULTS / "sRNA" / "tracks" / f"{name}__{size}__{strand}.bw"
                    assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_cluster_deg_results(self, results_exist):
        """Cluster DEG results for on_new_clusters and on_all_genes."""
        for mode in ["on_new_clusters", "on_all_genes"]:
            deg = RESULTS / "sRNA" / "clusters" / f"{RUN}__{mode}" / "DEG_WT_vs_dcr1.txt"
            assert _file_exists_nonempty(deg), f"Missing or empty: {deg}"

    def test_cluster_count_files(self, results_exist):
        """Cluster count files exist for both modes."""
        for mode in ["on_new_clusters", "on_all_genes"]:
            base = RESULTS / "sRNA" / "clusters" / f"{RUN}__{mode}"
            for fname in ["Counts.txt", "counts_for_edgeR.txt", "Results.txt"]:
                path = base / fname
                assert _file_exists_nonempty(path), f"Missing or empty: {path}"

    def test_srna_checkpoints(self, results_exist):
        """All sRNA checkpoints present."""
        for chkpt in SRNA_CHECKPOINTS:
            path = RESULTS / "sRNA" / "chkpts" / chkpt
            assert path.exists(), f"Missing checkpoint: {path}"


# ===========================================================================
# TestCombinedOutputs
# ===========================================================================

@pytest.mark.slow
@pytest.mark.integration
class TestCombinedOutputs:
    """Validate combined analysis outputs."""

    def test_heatmap_matrices(self, results_exist):
        """Heatmap matrices for tss/tes/regions exist."""
        for region_type in ["tss", "tes", "regions"]:
            gz = RESULTS / "combined" / "matrix" / f"final_matrix_{region_type}__most__{RUN}__all_genes.gz"
            assert _file_exists_nonempty(gz), f"Missing or empty: {gz}"

    def test_heatmap_matrices_have_one_region_group(self, results_exist):
        """A merged matrix carries exactly one region group.

        all_genes is a stranded target, so making_stranded_matrix_on_targetfile
        builds it in two halves and merging_matrix rbinds them. Everything
        downstream passes a single --regionsLabel, so the halves have to come
        back as one group; two groups means the plots below cannot be drawn.
        """
        for region_type in ["tss", "tes", "regions"]:
            gz = RESULTS / "combined" / "matrix" / f"final_matrix_{region_type}__most__{RUN}__all_genes.gz"
            with gzip.open(gz, "rt") as fh:
                header = json.loads(fh.readline().lstrip("@"))
            assert len(header["group_labels"]) == 1, (
                f"{gz.name} has {len(header['group_labels'])} region groups "
                f"({header['group_labels']}), expected 1"
            )

    def test_heatmap_pdfs(self, results_exist):
        """Heatmap PDFs for tss/tes/regions exist."""
        for region_type in ["tss", "tes", "regions"]:
            pdf = RESULTS / "combined" / "plots" / f"Heatmap__{region_type}__most__{RUN}__all_genes.pdf"
            assert _file_exists_nonempty(pdf), f"Missing or empty: {pdf}"

    def test_profile_pdfs(self, results_exist):
        """Profile PDFs for tss/tes/regions exist."""
        for region_type in ["tss", "tes", "regions"]:
            pdf = RESULTS / "combined" / "plots" / f"Profile__{region_type}__most__{RUN}__all_genes.pdf"
            assert _file_exists_nonempty(pdf), f"Missing or empty: {pdf}"

    def test_browser_track_pdf(self, results_exist):
        """Browser track PDF exists."""
        pdf = RESULTS / "combined" / "plots" / f"Browser_full_chromosomes__all__{RUN}.pdf"
        assert _file_exists_nonempty(pdf), f"Missing or empty: {pdf}"

    def test_pca_matrices(self, results_exist):
        """PCA matrices (.npz) exist."""
        for prefix in ["pca_matrix__all_chip", "pca_matrix__ChIP"]:
            npz = RESULTS / "combined" / "matrix" / f"{prefix}__{RUN}.npz"
            assert _file_exists_nonempty(npz), f"Missing or empty: {npz}"

    def test_bed_annotation_files(self, results_exist):
        """BED annotation files exist."""
        bedfiles = RESULTS / "combined" / "bedfiles"
        expected = [
            f"{GENOME}__all_genes.bed",
            f"{GENOME}__protein_coding_genes.bed",
            f"full_chromosomes__{GENOME}.bed",
            f"combined_peaks__ChIP__{RUN}.bed",
            f"annotated__combined_peaks__ChIP__{RUN}.bed",
        ]
        for fname in expected:
            path = bedfiles / fname
            assert _file_exists_nonempty(path), f"Missing or empty: {path}"


# ===========================================================================
# TestOutputIntegrity
# ===========================================================================

@pytest.mark.slow
@pytest.mark.integration
class TestOutputIntegrity:
    """Validate file format integrity of key outputs."""

    def test_bam_header_has_pombe_chromosomes(self, results_exist):
        """BAM header contains @SQ lines with pombe chromosomes (I, II, III)."""
        import subprocess
        bam = RESULTS / "ChIP" / "mapped" / f"final__{CHIP_PE_PEAK_SAMPLES[0]}.bam"
        result = subprocess.run(
            ["samtools", "view", "-H", str(bam)],
            capture_output=True, text=True, timeout=30
        )
        assert result.returncode == 0, f"samtools view -H failed: {result.stderr}"
        header = result.stdout
        assert "@SQ" in header, "No @SQ lines in BAM header"
        # S. pombe chromosomes are named I, II, III
        sq_lines = [l for l in header.splitlines() if l.startswith("@SQ")]
        chrom_names = [l.split("\t")[1].replace("SN:", "") for l in sq_lines]
        for expected in ["I", "II", "III"]:
            assert expected in chrom_names, \
                f"Chromosome {expected} not in BAM header: {chrom_names}"

    def test_bigwig_magic_bytes(self, results_exist):
        """Bigwig files have valid magic bytes (0x888FFC26 or 0x26FC8F88)."""
        bw = RESULTS / "ChIP" / "tracks" / f"FC__final__{CHIP_NON_INPUT_PER_REP[0]}.bw"
        magic = _read_bigwig_magic(bw)
        valid_magics = {0x888FFC26, 0x26FC8F88}
        assert magic in valid_magics, \
            f"Invalid bigwig magic: 0x{magic:08X}, expected one of {[hex(m) for m in valid_magics]}"

    def test_broadpeak_has_9_columns(self, results_exist):
        """BroadPeak files have 9 tab-separated columns."""
        peak = RESULTS / "ChIP" / "peaks" / f"peaks_pe__final__{CHIP_PE_PEAK_SAMPLES[0]}_peaks.broadPeak"
        with open(peak) as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                cols = line.split("\t")
                assert len(cols) == 9, \
                    f"BroadPeak line {i+1} has {len(cols)} columns, expected 9"
                break  # only check first data line

    def test_narrowpeak_has_10_columns(self, results_exist):
        """NarrowPeak files have 10 tab-separated columns."""
        peak = RESULTS / "ChIP" / "peaks" / f"peaks_se__final__{CHIP_SE_PEAK_SAMPLES[0]}_peaks.narrowPeak"
        with open(peak) as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                cols = line.split("\t")
                assert len(cols) == 10, \
                    f"NarrowPeak line {i+1} has {len(cols)} columns, expected 10"
                break  # only check first data line
