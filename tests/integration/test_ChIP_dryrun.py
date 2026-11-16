"""
Integration tests for ChIP-seq workflows using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for ChIP_broad and ChIP_narrow assay types without executing rules.

Test samples (from test_samples_ChIP.tsv):
  ChIP_broad (H3K9me2):
    WT_ChIP_H3K9me2_rep1        ChIP_broad   PE   genotype:WT       (ctrl: WT_ChIP_WCE_rep1)
    WT_ChIP_H3K9me2_rep2        ChIP_broad   PE   genotype:WT       (ctrl: WT_ChIP_WCE_rep2)
    mutant_ChIP_H3K9me2_rep1    ChIP_broad   PE   genotype:mutant   (ctrl: mutant_ChIP_WCE_rep1)
    mutant_ChIP_H3K9me2_rep2    ChIP_broad   PE   genotype:mutant   (ctrl: mutant_ChIP_WCE_rep1)
    WT_ChIP_WCE_rep1            ChIP_broad   PE   genotype:WT       (control)
    WT_ChIP_WCE_rep2            ChIP_broad   PE   genotype:WT       (control)
    mutant_ChIP_WCE_rep1        ChIP_broad   PE   genotype:mutant   (control)
    mutant_ChIP_WCE_rep2        ChIP_broad   PE   genotype:mutant   (control)
  ChIP_narrow (H3K4me3):
    WT_ChIP_H3K4me3_rep1        ChIP_narrow  SE   genotype:WT       (ctrl: WT_ChIP_Input_rep1)
    WT_ChIP_H3K4me3_rep2        ChIP_narrow  SE   genotype:WT       (ctrl: WT_ChIP_Input_rep1)
    mutant_ChIP_H3K4me3_rep1    ChIP_narrow  SE   genotype:mutant   (ctrl: WT_ChIP_Input_rep1)
    mutant_ChIP_H3K4me3_rep2    ChIP_narrow  SE   genotype:mutant   (ctrl: WT_ChIP_Input_rep1)
    WT_ChIP_Input_rep1          ChIP_narrow  SE   genotype:WT       (control)
    mutant_ChIP_Input_rep1      ChIP_narrow  SE   genotype:mutant   (control)
"""

import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun, run_snakemake_dag,
)


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_ChIP.yaml")


# ---------------------------------------------------------------------------
# Target path constants
# ---------------------------------------------------------------------------

# Per-replicate bigwig targets (FC__final__{sample_name}.bw)
BROAD_PE_BIGWIG = f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__WT_ChIP_H3K9me2_rep1.bw"
BROAD_PE_REP2_BIGWIG = f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__WT_ChIP_H3K9me2_rep2.bw"
NARROW_SE_BIGWIG = f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__WT_ChIP_H3K4me3_rep1.bw"
MUTANT_BROAD_BIGWIG = f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__mutant_ChIP_H3K9me2_rep1.bw"

# Per-replicate peak targets
BROAD_PE_PEAKS = f"{_OUTPUT_DIR}/ChIP/peaks/peaks_pe__final__WT_ChIP_H3K9me2_rep1_peaks.broadPeak"
NARROW_SE_PEAKS = f"{_OUTPUT_DIR}/ChIP/peaks/peaks_se__final__WT_ChIP_H3K4me3_rep1_peaks.narrowPeak"

# Analysis-level names (Assay__levels_label__IP_target__Genome)
BROAD_WT_ANALYSIS = "ChIP_broad__WT__H3K9me2__test_genome"
BROAD_MUTANT_ANALYSIS = "ChIP_broad__mutant__H3K9me2__test_genome"
NARROW_WT_ANALYSIS = "ChIP_narrow__WT__H3K4me3__test_genome"

# Analysis-level targets
BROAD_MERGED_BIGWIG = f"{_OUTPUT_DIR}/ChIP/tracks/FC__merged__{BROAD_WT_ANALYSIS}.bw"
BROAD_SELECTED_PEAKS = f"{_OUTPUT_DIR}/ChIP/peaks/selected_peaks__{BROAD_WT_ANALYSIS}.bedPeak"
NARROW_SELECTED_PEAKS = f"{_OUTPUT_DIR}/ChIP/peaks/selected_peaks__{NARROW_WT_ANALYSIS}.bedPeak"

# Differential peaks (MAnorm) target
DIFF_PEAKS_TARGET = f"{_OUTPUT_DIR}/ChIP/peaks/{BROAD_WT_ANALYSIS}_vs_{BROAD_MUTANT_ANALYSIS}/{BROAD_WT_ANALYSIS}_vs_{BROAD_MUTANT_ANALYSIS}_all_MAvalues.xls"


@pytest.fixture(scope="module")
def test_options(repo_root):
    """Get the path to the test options file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_ChIP.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


# ===========================================================================
# TestBasic
# ===========================================================================

class TestBasic:
    """Basic sanity checks for the ChIP test configuration."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        """Test that the test options file exists."""
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_ChIP.tsv"
        assert sample_file.exists()


# ===========================================================================
# TestChIPBroadWorkflow
# ===========================================================================

class TestChIPBroadWorkflow:
    """Test ChIP_broad (H3K9me2) PE workflow dry-run."""

    def test_broad_pe_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ChIP_broad PE bigwig target."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BROAD_PE_BIGWIG)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_broad_pe_includes_filter_bam_pe(self, snakemake_available, repo_root, test_options):
        """Test that ChIP_broad PE workflow includes filter_bam_pe rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BROAD_PE_BIGWIG, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "filter_bam_pe" in output, "ChIP_broad PE should use filter_bam_pe rule"

    def test_broad_pe_peaks_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ChIP_broad PE peak target."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BROAD_PE_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# TestChIPNarrowWorkflow
# ===========================================================================

class TestChIPNarrowWorkflow:
    """Test ChIP_narrow (H3K4me3) SE workflow dry-run."""

    def test_narrow_se_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ChIP_narrow SE bigwig target."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, NARROW_SE_BIGWIG)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_narrow_se_includes_filter_bam_se(self, snakemake_available, repo_root, test_options):
        """Test that ChIP_narrow SE workflow includes filter_bam_se rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, NARROW_SE_BIGWIG, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "filter_bam_se" in output, "ChIP_narrow SE should use filter_bam_se rule"

    def test_narrow_se_peaks_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ChIP_narrow SE peak target."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, NARROW_SE_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# TestReplicateMerging
# ===========================================================================

class TestReplicateMerging:
    """Test replicate merging rules for assays with multiple replicates."""

    def test_broad_two_reps_trigger_merging(self, snakemake_available, repo_root, test_options):
        """WT H3K9me2 has 2 reps -- merging_bam_replicates should appear."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BROAD_MERGED_BIGWIG, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "merging_bam_replicates" in output, \
            "ChIP_broad with 2 reps should trigger merging_bam_replicates rule"

    def test_selected_peaks_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that selected_peaks target for broad WT succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BROAD_SELECTED_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_narrow_selected_peaks_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that selected_peaks target for narrow WT succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, NARROW_SELECTED_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# TestDifferentialPeaks
# ===========================================================================

class TestDifferentialPeaks:
    """Test MAnorm pairwise differential peak analysis (WT vs mutant H3K9me2)."""

    def test_diff_peaks_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that MAnorm dry-run succeeds for WT vs mutant H3K9me2."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DIFF_PEAKS_TARGET)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_diff_peaks_includes_manorm(self, snakemake_available, repo_root, test_options):
        """Test that differential peaks workflow includes perform_pairwise_diff_peaks."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DIFF_PEAKS_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "perform_pairwise_diff_peaks" in output, \
            "Differential peaks should include perform_pairwise_diff_peaks rule"


# ===========================================================================
# TestWildcardResolution
# ===========================================================================

class TestWildcardResolution:
    """Test wildcard resolution for various ChIP sample and target combinations."""

    def test_broad_pe_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for broad PE samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__WT_ChIP_H3K9me2_rep1.bw",
            f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__WT_ChIP_H3K9me2_rep2.bw",
            f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__mutant_ChIP_H3K9me2_rep1.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_narrow_se_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for narrow SE samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__WT_ChIP_H3K4me3_rep1.bw",
            f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__WT_ChIP_H3K4me3_rep2.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_peak_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that peak file wildcards resolve for both broad and narrow."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/ChIP/peaks/peaks_pe__final__WT_ChIP_H3K9me2_rep1_peaks.broadPeak",
            f"{_OUTPUT_DIR}/ChIP/peaks/peaks_se__final__WT_ChIP_H3K4me3_rep1_peaks.narrowPeak",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_analysis_level_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that analysis-level targets resolve."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            BROAD_MERGED_BIGWIG,
            BROAD_SELECTED_PEAKS,
            NARROW_SELECTED_PEAKS,
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"


# ===========================================================================
# TestDAGStructure
# ===========================================================================

class TestDAGStructure:
    """Test DAG generation and structure for ChIP-seq rules."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, BROAD_PE_BIGWIG)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_chip_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG for a ChIP target contains ChIP-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, BROAD_PE_BIGWIG)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        for rule in ("filter_bam_pe", "dispatch_final_bam", "make_bigwig_chip"):
            assert rule in dag_output, f"Rule '{rule}' not found in ChIP DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_options):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, BROAD_PE_BIGWIG)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert "->" in result.stdout, "DAG should contain rule dependencies"


# ===========================================================================
# TestErrorHandling
# ===========================================================================

class TestErrorHandling:
    """Test error handling for invalid configurations."""

    def test_missing_sample_target(self, snakemake_available, repo_root, test_options):
        """Test that requesting a target for a non-existent sample fails gracefully."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/ChIP/tracks/FC__final__nonexistent_sample.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for non-existent sample"
