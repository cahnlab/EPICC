"""
Integration tests for S. pombe multi-assay workflow using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for all four assay types without actually executing the rules.

Test samples (from test_samples_pombe.tsv):
  ChIP_broad (7): WT/dcr1 x H3K9me2/me3 + 1 shared WCE, PE
  ChIP_narrow (3): veg H3K4me3 x 2 reps + Input, SE
  RNAseq (4): WT/dcr1 x 2 reps, SE
  sRNA (3): WT x 2 reps + dcr1 x 1 rep, SE

Output paths use Sample_ID directly for per-replicate files.
Analysis-level names use build_analysis_name():
  {Assay}__{levels_label}__{IP_target}__{Genome} (empty parts omitted)
"""

import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun, run_snakemake_dag,
)


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_pombe.yaml")


# ---------------------------------------------------------------------------
# Target path constants (derived from test_samples_pombe.tsv)
# ---------------------------------------------------------------------------

# Per-replicate targets
CHIP_BROAD_TARGET = f"{_OUTPUT_DIR}/ChIP/tracks/coverage__WT_ChIP_H3K9me2_rep1.bw"
CHIP_BROAD_TARGET_REP2 = f"{_OUTPUT_DIR}/ChIP/tracks/coverage__WT_ChIP_H3K9me2_rep2.bw"
CHIP_NARROW_TARGET = f"{_OUTPUT_DIR}/ChIP/tracks/coverage__WT_ChIP_H3K4me3_rep1.bw"
RNA_TARGET = f"{_OUTPUT_DIR}/RNA/mapped/final__WT_RNA_rep1.bam"
SRNA_TARGET = f"{_OUTPUT_DIR}/sRNA/mapped/WT_sRNA_rep1/Results.txt"

# Control sample targets
INPUT_BROAD_TARGET = f"{_OUTPUT_DIR}/ChIP/tracks/coverage__WT_ChIP_WCE_rep1.bw"
INPUT_NARROW_TARGET = f"{_OUTPUT_DIR}/ChIP/tracks/coverage__WT_ChIP_Input_rep1.bw"

# Analysis-level names (Assay__levels_label__IP_target__Genome)
CHIP_BROAD_ANALYSIS = "ChIP_broad__WT__H3K9me2__Spombe"
CHIP_NARROW_ANALYSIS = "ChIP_narrow__WT__H3K4me3__Spombe"
RNA_ANALYSIS_WT = "RNAseq__WT__Spombe"
SRNA_ANALYSIS_WT = "sRNA__WT__Spombe"

# PCA matrix over the ChIP BAM tracks (issue #50)
CHIP_PCA_TARGET = f"{_OUTPUT_DIR}/combined/matrix/pca_matrix__ChIP__test_pombe__Spombe.npz"

# Env checkpoint targets
CHIP_CHECKPOINT = f"{_OUTPUT_DIR}/ChIP/chkpts/ChIP_analysis__test_pombe__Spombe.done"
RNA_CHECKPOINT = f"{_OUTPUT_DIR}/RNA/chkpts/RNA_analysis__test_pombe__Spombe.done"
SRNA_CHECKPOINT = f"{_OUTPUT_DIR}/sRNA/chkpts/sRNA_analysis__test_pombe__Spombe.done"
FINAL_CHECKPOINT = f"{_OUTPUT_DIR}/combined/chkpts/final_analysis__test_pombe.done"


@pytest.fixture(scope="module")
def test_options(repo_root):
    """Get the path to the test options file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_pombe.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


class TestPombeDryRunBasic:
    """Basic setup verification tests."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        """Test that the test options file exists."""
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_pombe.tsv"
        assert sample_file.exists()


class TestChIPBroadWorkflow:
    """Test ChIP-seq broad-peak workflow dry-run (PE samples)."""

    def test_chip_broad_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ChIP broad coverage bigwig."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_BROAD_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_chip_broad_includes_expected_rules(self, snakemake_available, repo_root, test_options):
        """Test that ChIP broad workflow includes PE mapping and broad peak rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_BROAD_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "filter_bam_pe",
            "make_coverage_chip",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_chip_broad_excludes_se_mapping(self, snakemake_available, repo_root, test_options):
        """Test that PE ChIP broad samples do not trigger SE mapping rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_BROAD_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "filter_bam_se" not in output, \
            "SE mapping rule should not appear for PE ChIP broad sample"

    def test_chip_broad_peak_calling(self, snakemake_available, repo_root, test_options):
        """Test that ChIP broad uses PE peak calling with broad mode."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "calling_peaks_pe" in output, \
            "PE peak calling rule should be in ChIP broad workflow"
        assert "--broad" in output, \
            "ChIP broad should use --broad flag in MACS2"


class TestChIPNarrowWorkflow:
    """Test ChIP-seq narrow-peak workflow dry-run (SE samples)."""

    def test_chip_narrow_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ChIP narrow coverage bigwig."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_NARROW_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_chip_narrow_includes_expected_rules(self, snakemake_available, repo_root, test_options):
        """Test that ChIP narrow workflow includes SE mapping rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_NARROW_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "filter_bam_se",
            "make_coverage_chip",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_chip_narrow_excludes_pe_mapping(self, snakemake_available, repo_root, test_options):
        """Test that SE ChIP narrow samples do not trigger PE mapping rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_NARROW_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "filter_bam_pe" not in output, \
            "PE mapping rule should not appear for SE ChIP narrow sample"


class TestRNAseqWorkflow:
    """Test RNA-seq workflow dry-run (SE samples)."""

    def test_rnaseq_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for RNA-seq final BAM."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, RNA_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_rnaseq_includes_expected_rules(self, snakemake_available, repo_root, test_options):
        """Test that RNA-seq workflow includes SE STAR mapping and filtering."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, RNA_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "STAR_map_se",
            "filter_rna_se",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_rnaseq_deg_target(self, snakemake_available, repo_root, test_options):
        """Test that DEG calling dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, RNA_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "call_all_DEGs" in output, \
            "RNA checkpoint should include DEG calling rule"


class TestSmallRNAWorkflow:
    """Test small RNA-seq workflow dry-run (SE samples)."""

    def test_srna_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ShortStack results."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, SRNA_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_srna_includes_expected_rules(self, snakemake_available, repo_root, test_options):
        """Test that sRNA workflow includes ShortStack and cluster rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, SRNA_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "shortstack_map",
            "merging_srna_replicates",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_srna_structural_rna_depletion(self, snakemake_available, repo_root, test_options):
        """Test that structural RNA depletion rules are present."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, SRNA_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "filter_structural_rna" in output, \
            "Structural RNA depletion should be present (structural_rna_depletion: true)"


class TestControlLinking:
    """Test that ChIP samples correctly resolve their Input controls."""

    def test_chip_broad_resolves_control(self, snakemake_available, repo_root, test_options):
        """Test that ChIP broad workflow references WT Input control in peak calling."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "WT_ChIP_WCE_rep1" in output, \
            "ChIP broad workflow should reference WT_ChIP_WCE_rep1 as control"

    def test_chip_narrow_resolves_control(self, snakemake_available, repo_root, test_options):
        """Test that ChIP narrow workflow references WT Input control in peak calling."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "WT_ChIP_Input_rep1" in output, \
            "ChIP workflow should reference WT_ChIP_Input_rep1 as control"

    def test_input_control_processed(self, snakemake_available, repo_root, test_options):
        """Test that Input control samples can be processed independently."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, INPUT_BROAD_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


class TestDAGStructure:
    """Test DAG generation and structure."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DAG can be generated for final checkpoint."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_options, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_all_assay_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG contains rules from all four assay types."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_options, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        dag_output = result.stdout

        # Rules from each assay type
        assay_rules = {
            "ChIP": "filter_bam",
            "RNA": "STAR_map",
            "sRNA": "shortstack_map",
        }

        for assay, rule in assay_rules.items():
            assert rule in dag_output, f"Rule '{rule}' for {assay} not found in DAG"

    def test_dag_has_dependencies(self, snakemake_available, repo_root, test_options):
        """Test that DAG has rule dependencies (edges)."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_options, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        assert "->" in result.stdout, "DAG should contain rule dependencies"


class TestReplicateHandling:
    """Test handling of multiple replicates."""

    def test_chip_broad_replicates_resolve(self, snakemake_available, repo_root, test_options):
        """Test that both ChIP broad H3K9me2 replicates can be resolved."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        for target in [CHIP_BROAD_TARGET, CHIP_BROAD_TARGET_REP2]:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Replicate failed for {target}: {result.stderr}"

    def test_idr_analysis_present(self, snakemake_available, repo_root, test_options):
        """Test that IDR analysis rule is present for replicated ChIP samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "idr_analysis_replicates" in output, \
            "IDR analysis should be present for replicated ChIP samples"

    def test_rna_replicates_resolve(self, snakemake_available, repo_root, test_options):
        """Test that multiple RNA-seq replicates can be resolved."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        targets = [
            f"{_OUTPUT_DIR}/RNA/mapped/final__WT_RNA_rep1.bam",
            f"{_OUTPUT_DIR}/RNA/mapped/final__WT_RNA_rep2.bam",
        ]

        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Replicate failed for {target}: {result.stderr}"


class TestEnvCheckpoints:
    """Test environment-level checkpoint targets (validates full DAG resolution)."""

    def test_chip_checkpoint(self, snakemake_available, repo_root, test_options):
        """Test that ChIP analysis checkpoint dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, CHIP_CHECKPOINT)

        assert result.returncode == 0, f"ChIP checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_rna_checkpoint(self, snakemake_available, repo_root, test_options):
        """Test that RNA analysis checkpoint dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, RNA_CHECKPOINT)

        assert result.returncode == 0, f"RNA checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_srna_checkpoint(self, snakemake_available, repo_root, test_options):
        """Test that sRNA analysis checkpoint dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, SRNA_CHECKPOINT)

        assert result.returncode == 0, f"sRNA checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_final_checkpoint(self, snakemake_available, repo_root, test_options):
        """Test that final combined checkpoint dry-run succeeds.

        This is the most comprehensive test — it forces Snakemake to resolve
        the entire multi-assay DAG (~257 steps) in dry-run mode.
        """
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"Final checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


class TestMapOnly:
    """Test intermediate target rules."""

    def test_map_only_target(self, snakemake_available, repo_root, test_options):
        """Test that map_only target dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, "map_only")

        assert result.returncode == 0, f"map_only failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_coverage_chip_target(self, snakemake_available, repo_root, test_options):
        """Test that coverage_chip target dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, "coverage_chip")

        assert result.returncode == 0, f"coverage_chip failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


def _pca_rule_inputs(debug_dag_output):
    """Extract the input file list of the summarize_tracks_pca job.

    `--debug-dag` prints each resolved job as a `rule <name>:` block whose
    `input:` line is a comma-separated file list.
    """
    lines = debug_dag_output.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("rule summarize_tracks_pca:"):
            for follow in lines[i + 1:i + 6]:
                stripped = follow.strip()
                if stripped.startswith("input:"):
                    return [f.strip() for f in
                            stripped[len("input:"):].split(",") if f.strip()]
    return []


class TestPCATrackIndexes:
    """PCA over BAM tracks must declare the .bai indexes (issue #50).

    multiBamSummary reads each BAM's index, but the shell command never names
    it. Without a declared input edge Snakemake can neither rebuild a missing
    index nor keep a temp() one alive for this job, which surfaced as
    intermittent "missing an index" failures in summarize_tracks_pca.
    """

    def test_chip_pca_declares_bam_indexes(self, snakemake_available, repo_root, test_options):
        """Every BAM input of summarize_tracks_pca has its .bai declared too."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(
            repo_root, test_options, CHIP_PCA_TARGET, ["--debug-dag"]
        )
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        inputs = _pca_rule_inputs(result.stdout + result.stderr)
        assert inputs, "summarize_tracks_pca job not found in the DAG"

        bams = [f for f in inputs if f.endswith(".bam")]
        bais = {f for f in inputs if f.endswith(".bam.bai")}
        assert bams, "expected BAM tracks for the ChIP PCA"
        missing = [b for b in bams if b + ".bai" not in bais]
        assert not missing, f"BAM inputs without a declared index: {missing}"

    def test_chip_pca_inputs_are_only_bams_and_indexes(self, snakemake_available, repo_root, test_options):
        """The index list is scoped to BAM envs, not appended blindly.

        The mC contexts summarize bigwigs and must get no indexes; pombe has no
        mC samples, so assert the ChIP job's inputs are exclusively BAM/BAI.
        """
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(
            repo_root, test_options, CHIP_PCA_TARGET, ["--debug-dag"]
        )
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        inputs = _pca_rule_inputs(result.stdout + result.stderr)
        assert inputs, "summarize_tracks_pca job not found in the DAG"
        assert all(f.endswith((".bam", ".bam.bai")) for f in inputs), \
            f"unexpected non-BAM inputs for a ChIP PCA: {inputs}"
