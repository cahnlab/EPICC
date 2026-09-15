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

import pytest
from pathlib import Path

from tests.integration.conftest import REPO_ROOT, load_output_dir

# ---------------------------------------------------------------------------
# Constants: repo root and results directory
# ---------------------------------------------------------------------------
RESULTS = REPO_ROOT / load_output_dir("test_options_pombe.yaml")

# ---------------------------------------------------------------------------
# Expected file lists (verified against actual pombe run 2026-02-26)
# ---------------------------------------------------------------------------

# ChIP per-replicate BAMs (10 samples)
CHIP_PER_REP_BAMS = [
    "WT_H3K9me2_rep1",
    "WT_H3K9me2_rep2",
    "WT_H3K9me3_rep1",
    "WT_H3K9me3_rep2",
    "dcr1_H3K9me2_rep1",
    "dcr1_H3K9me3_rep1",
    "WT_WCE_rep1",
    "WT_H3K4me3_rep1",
    "WT_H3K4me3_rep2",
    "WT_Input_rep1",
]

# ChIP analysis-level merged BAMs (3 analysis groups + 2 merged controls)
CHIP_MERGED_ANALYSIS_BAMS = [
    "ChIP_broad__WT__H3K9me2__Spombe",
    "ChIP_broad__WT__H3K9me3__Spombe",
    "ChIP_narrow__WT__H3K4me3__Spombe",
]
CHIP_MERGED_CONTROL_BAMS = [
    "WT_WCE_rep1",
    "WT_Input_rep1",
]

# ChIP non-Input/WCE samples (for FC bigwigs)
CHIP_NON_INPUT_PER_REP = [
    "WT_H3K9me2_rep1",
    "WT_H3K9me2_rep2",
    "WT_H3K9me3_rep1",
    "WT_H3K9me3_rep2",
    "dcr1_H3K9me2_rep1",
    "dcr1_H3K9me3_rep1",
    "WT_H3K4me3_rep1",
    "WT_H3K4me3_rep2",
]

# ChIP PE samples that should have broadPeak files
CHIP_PE_PEAK_SAMPLES = [
    "WT_H3K9me2_rep1",
    "WT_H3K9me2_rep2",
    "WT_H3K9me3_rep1",
    "WT_H3K9me3_rep2",
    "dcr1_H3K9me2_rep1",
    "dcr1_H3K9me3_rep1",
]

# ChIP SE samples that should have narrowPeak files
CHIP_SE_PEAK_SAMPLES = [
    "WT_H3K4me3_rep1",
    "WT_H3K4me3_rep2",
]

# RNA per-replicate samples (4)
RNA_PER_REP = [
    "WT_RNA_rep1",
    "WT_RNA_rep2",
    "dcr1_RNA_rep1",
    "dcr1_RNA_rep2",
]

# RNA merged BAMs (2)
RNA_MERGED = [
    "RNAseq__WT__Spombe",
    "RNAseq__dcr1__Spombe",
]

# sRNA per-replicate samples (3)
SRNA_PER_REP = [
    "WT_sRNA_rep1",
    "WT_sRNA_rep2",
    "dcr1_sRNA_rep1",
]

# sRNA size classes
SRNA_SIZES = ["21nt", "22nt", "23nt", "24nt"]
SRNA_STRANDS = ["plus", "minus"]

# Checkpoint files
CHIP_CHECKPOINTS = [
    "ChIP_analysis__test_pombe__Spombe.done",
    "idr__ChIP_broad__WT__H3K9me2__Spombe.done",
    "idr__ChIP_broad__WT__H3K9me3__Spombe.done",
    "idr__ChIP_narrow__WT__H3K4me3__Spombe.done",
] + [f"map_ChIP__{s}.done" for s in CHIP_PER_REP_BAMS]

RNA_CHECKPOINTS = [
    "calling_DEGs__test_pombe__Spombe.done",
    "RNA_analysis__test_pombe__Spombe.done",
] + [f"map_RNA__{s}.done" for s in RNA_PER_REP]

SRNA_CHECKPOINTS = [
    "sRNA_analysis__test_pombe__Spombe.done",
    "calling_differential_sRNA_clusters__test_pombe__Spombe__on_all_genes.done",
    "calling_differential_sRNA_clusters__test_pombe__Spombe__on_new_clusters.done",
] + [f"map_sRNA__{s}.done" for s in SRNA_PER_REP]

COMBINED_CHECKPOINTS = [
    "combined_analysis__test_pombe__Spombe.done",
    "final_analysis__test_pombe.done",
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
        """All 11 per-replicate BAMs and BAIs exist and are non-empty."""
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
        """FC bigwigs for 8 non-Input per-rep samples exist."""
        for sample in CHIP_NON_INPUT_PER_REP:
            bw = RESULTS / "ChIP" / "tracks" / f"FC__final__{sample}.bw"
            assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_fc_bigwigs_merged(self, results_exist):
        """Merged FC bigwigs for 3 analysis groups exist."""
        for name in CHIP_MERGED_ANALYSIS_BAMS:
            bw = RESULTS / "ChIP" / "tracks" / f"FC__merged__{name}.bw"
            assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_broad_peaks(self, results_exist):
        """BroadPeak files for 6 PE samples exist."""
        for sample in CHIP_PE_PEAK_SAMPLES:
            peak = RESULTS / "ChIP" / "peaks" / f"peaks_pe__final__{sample}_peaks.broadPeak"
            assert _file_exists_nonempty(peak), f"Missing or empty: {peak}"

    def test_narrow_peaks(self, results_exist):
        """NarrowPeak files for 2 SE samples exist."""
        for sample in CHIP_SE_PEAK_SAMPLES:
            peak = RESULTS / "ChIP" / "peaks" / f"peaks_se__final__{sample}_peaks.narrowPeak"
            assert _file_exists_nonempty(peak), f"Missing or empty: {peak}"

    def test_idr_files(self, results_exist):
        """IDR output files exist (may be empty when IDR fails due to few peaks)."""
        for name in ["ChIP_broad__WT__H3K9me2__Spombe",
                      "ChIP_broad__WT__H3K9me3__Spombe"]:
            idr = RESULTS / "ChIP" / "peaks" / f"idr_peaks__{name}.broadPeak"
            assert idr.exists(), f"Missing: {idr}"
        idr_narrow = RESULTS / "ChIP" / "peaks" / "idr_peaks__ChIP_narrow__WT__H3K4me3__Spombe.narrowPeak"
        assert idr_narrow.exists(), f"Missing: {idr_narrow}"

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
        """4 per-replicate BAMs + BAIs exist."""
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
        rpkm = RESULTS / "RNA" / "DEG" / "genes_rpkm__test_pombe__Spombe.txt"
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
        deg = RESULTS / "RNA" / "DEG" / "DEG_test_pombe__Spombe__WT___vs_dcr1__.txt"
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
        """ShortStack Results.txt for 3 samples."""
        for sample in SRNA_PER_REP:
            results_txt = RESULTS / "sRNA" / "mapped" / sample / "Results.txt"
            assert _file_exists_nonempty(results_txt), f"Missing or empty: {results_txt}"

    def test_size_class_bigwigs_per_rep(self, results_exist):
        """Size-class bigwigs (4 sizes x 2 strands x 3 reps)."""
        for sample in SRNA_PER_REP:
            for size in SRNA_SIZES:
                for strand in SRNA_STRANDS:
                    bw = RESULTS / "sRNA" / "tracks" / f"{sample}__{size}__{strand}.bw"
                    assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_merged_wt_bigwigs(self, results_exist):
        """Merged WT bigwigs (4 sizes x 2 strands)."""
        for size in SRNA_SIZES:
            for strand in SRNA_STRANDS:
                bw = RESULTS / "sRNA" / "tracks" / f"sRNA__WT__Spombe__{size}__{strand}.bw"
                assert _file_exists_nonempty(bw), f"Missing or empty: {bw}"

    def test_cluster_deg_results(self, results_exist):
        """Cluster DEG results for on_new_clusters and on_all_genes."""
        for mode in ["on_new_clusters", "on_all_genes"]:
            deg = RESULTS / "sRNA" / "clusters" / f"test_pombe__Spombe__{mode}" / "DEG_WT_vs_dcr1.txt"
            assert _file_exists_nonempty(deg), f"Missing or empty: {deg}"

    def test_cluster_count_files(self, results_exist):
        """Cluster count files exist for both modes."""
        for mode in ["on_new_clusters", "on_all_genes"]:
            base = RESULTS / "sRNA" / "clusters" / f"test_pombe__Spombe__{mode}"
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
            gz = RESULTS / "combined" / "matrix" / f"final_matrix_{region_type}__most__test_pombe__Spombe__all_genes.gz"
            assert _file_exists_nonempty(gz), f"Missing or empty: {gz}"

    def test_heatmap_matrices_have_one_region_group(self, results_exist):
        """A merged matrix carries exactly one region group.

        all_genes is a stranded target, so making_stranded_matrix_on_targetfile
        builds it in two halves and merging_matrix rbinds them. Everything
        downstream passes a single --regionsLabel, so the halves have to come
        back as one group; two groups means the plots below cannot be drawn.
        """
        for region_type in ["tss", "tes", "regions"]:
            gz = RESULTS / "combined" / "matrix" / f"final_matrix_{region_type}__most__test_pombe__Spombe__all_genes.gz"
            with gzip.open(gz, "rt") as fh:
                header = json.loads(fh.readline().lstrip("@"))
            assert len(header["group_labels"]) == 1, (
                f"{gz.name} has {len(header['group_labels'])} region groups "
                f"({header['group_labels']}), expected 1"
            )

    def test_heatmap_pdfs(self, results_exist):
        """Heatmap PDFs for tss/tes/regions exist."""
        for region_type in ["tss", "tes", "regions"]:
            pdf = RESULTS / "combined" / "plots" / f"Heatmap__{region_type}__most__test_pombe__Spombe__all_genes.pdf"
            assert _file_exists_nonempty(pdf), f"Missing or empty: {pdf}"

    def test_profile_pdfs(self, results_exist):
        """Profile PDFs for tss/tes/regions exist."""
        for region_type in ["tss", "tes", "regions"]:
            pdf = RESULTS / "combined" / "plots" / f"Profile__{region_type}__most__test_pombe__Spombe__all_genes.pdf"
            assert _file_exists_nonempty(pdf), f"Missing or empty: {pdf}"

    def test_browser_track_pdf(self, results_exist):
        """Browser track PDF exists."""
        pdf = RESULTS / "combined" / "plots" / "Browser_full_chromosomes__all__test_pombe__Spombe.pdf"
        assert _file_exists_nonempty(pdf), f"Missing or empty: {pdf}"

    def test_pca_matrices(self, results_exist):
        """PCA matrices (.npz) exist."""
        for prefix in ["pca_matrix__all_chip", "pca_matrix__ChIP"]:
            npz = RESULTS / "combined" / "matrix" / f"{prefix}__test_pombe__Spombe.npz"
            assert _file_exists_nonempty(npz), f"Missing or empty: {npz}"

    def test_bed_annotation_files(self, results_exist):
        """BED annotation files exist."""
        bedfiles = RESULTS / "combined" / "bedfiles"
        expected = [
            "Spombe__all_genes.bed",
            "Spombe__protein_coding_genes.bed",
            "full_chromosomes__Spombe.bed",
            "combined_peaks__ChIP__test_pombe__Spombe.bed",
            "annotated__combined_peaks__ChIP__test_pombe__Spombe.bed",
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
        bam = RESULTS / "ChIP" / "mapped" / "final__WT_H3K9me2_rep1.bam"
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
        bw = RESULTS / "ChIP" / "tracks" / "FC__final__WT_H3K9me2_rep1.bw"
        magic = _read_bigwig_magic(bw)
        valid_magics = {0x888FFC26, 0x26FC8F88}
        assert magic in valid_magics, \
            f"Invalid bigwig magic: 0x{magic:08X}, expected one of {[hex(m) for m in valid_magics]}"

    def test_broadpeak_has_9_columns(self, results_exist):
        """BroadPeak files have 9 tab-separated columns."""
        peak = RESULTS / "ChIP" / "peaks" / "peaks_pe__final__WT_H3K9me2_rep1_peaks.broadPeak"
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
        peak = RESULTS / "ChIP" / "peaks" / "peaks_se__final__WT_H3K4me3_rep1_peaks.narrowPeak"
        with open(peak) as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                cols = line.split("\t")
                assert len(cols) == 10, \
                    f"NarrowPeak line {i+1} has {len(cols)} columns, expected 10"
                break  # only check first data line
