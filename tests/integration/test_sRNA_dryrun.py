"""
Integration tests for small RNA (sRNA) workflow using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for the sRNA assay type without actually executing the rules.

Test samples (from test_samples_sRNA.tsv):
  WT_sRNA_rep1        sRNA  SE  genotype:WT
  WT_sRNA_rep2        sRNA  SE  genotype:WT
  mutant_sRNA_rep1    sRNA  SE  genotype:mutant
  mutant_sRNA_rep2    sRNA  SE  genotype:mutant
  mutant2_sRNA_rep1   sRNA  SE  genotype:mutant2
"""

import os
import re
import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun, run_snakemake_dag,
)


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_sRNA.yaml")


# ---------------------------------------------------------------------------
# Target path constants
# ---------------------------------------------------------------------------

# Per-replicate targets are genome-qualified (mapped_name), as the pipeline
# builds them. A bare path is still constructible from the wildcards but is
# not what any generator requests.
_G = "test_genome"

# Per-replicate sized bigwig targets
WT_REP1_21NT_PLUS = f"{_OUTPUT_DIR}/sRNA/tracks/WT_sRNA_rep1__{_G}__21nt__plus.bw"
WT_REP1_24NT_PLUS = f"{_OUTPUT_DIR}/sRNA/tracks/WT_sRNA_rep1__{_G}__24nt__plus.bw"
WT_REP1_21NT_MINUS = f"{_OUTPUT_DIR}/sRNA/tracks/WT_sRNA_rep1__{_G}__21nt__minus.bw"

# Size stats target
WT_REP1_SIZE_STATS = f"{_OUTPUT_DIR}/sRNA/reports/sizes_stats__WT_sRNA_rep1__{_G}.txt"

# Cluster BED file target
WT_REP1_CLUSTERS = f"{_OUTPUT_DIR}/sRNA/mapped/WT_sRNA_rep1__{_G}/clusters.bed"

# Analysis-level names (Assay__levels_label__Genome)
SRNA_WT_ANALYSIS = "sRNA__WT__test_genome"
SRNA_MUTANT_ANALYSIS = "sRNA__mutant__test_genome"

# Merged-replicate bigwig targets (analysis-level)
WT_MERGED_21NT_PLUS = f"{_OUTPUT_DIR}/sRNA/tracks/{SRNA_WT_ANALYSIS}__21nt__plus.bw"
MUTANT_MERGED_21NT_PLUS = f"{_OUTPUT_DIR}/sRNA/tracks/{SRNA_MUTANT_ANALYSIS}__21nt__plus.bw"


@pytest.fixture(scope="module")
def test_options(repo_root):
    """Get the path to the test options file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_sRNA.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


# ===========================================================================
# TestBasic
# ===========================================================================

class TestBasic:
    """Basic sanity checks for the sRNA test configuration."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        """Test that the test options file exists."""
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_sRNA.tsv"
        assert sample_file.exists()


# ===========================================================================
# TestSRNAWorkflow
# ===========================================================================

class TestSRNAWorkflow:
    """Test sRNA workflow dry-run."""

    def test_srna_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for sRNA sized bigwig target."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_PLUS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_srna_includes_shortstack_map(self, snakemake_available, repo_root, test_options):
        """Test that sRNA workflow includes shortstack_map rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "shortstack_map" in output, "sRNA workflow should use shortstack_map rule"

    def test_srna_includes_expected_rules(self, snakemake_available, repo_root, test_options):
        """Test that sRNA workflow includes key expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        expected_rules = [
            "shortstack_map",
            "filter_size_srna_sample",
            "make_srna_stranded_bigwigs",
        ]
        for rule in expected_rules:
            assert rule in output, f"Expected sRNA rule '{rule}' not found in dry-run output"


# ===========================================================================
# TestStructuralRNADepletion
# ===========================================================================

class TestStructuralRNADepletion:
    """Test structural RNA depletion is in the DAG when enabled."""

    def test_filter_structural_rna_in_dag(self, snakemake_available, repo_root, test_options):
        """Test that filter_structural_rna rule is in the DAG."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "filter_structural_rna" in output, \
            "Structural RNA depletion should be in the DAG when structural_rna_depletion is true"

    def test_bt2_structural_index_in_dag(self, snakemake_available, repo_root, test_options):
        """Test that bowtie2 structural RNA index building is in the DAG."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "make_bt2_indices_for_structural_RNAs" in output, \
            "Structural RNA bt2 index should be in the DAG"


# ===========================================================================
# TestSizeFiltering
# ===========================================================================

class TestSizeFiltering:
    """Test size-filtered bigwig generation."""

    def test_21nt_plus_bigwig_resolves(self, snakemake_available, repo_root, test_options):
        """Test that 21nt plus-strand bigwig target resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_PLUS)
        assert result.returncode == 0, f"Dry-run failed for 21nt plus: {result.stderr}"

    def test_24nt_plus_bigwig_resolves(self, snakemake_available, repo_root, test_options):
        """Test that 24nt plus-strand bigwig target resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_24NT_PLUS)
        assert result.returncode == 0, f"Dry-run failed for 24nt plus: {result.stderr}"

    def test_21nt_minus_bigwig_resolves(self, snakemake_available, repo_root, test_options):
        """Test that 21nt minus-strand bigwig target resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_MINUS)
        assert result.returncode == 0, f"Dry-run failed for 21nt minus: {result.stderr}"

    def test_filter_size_rule_in_dag(self, snakemake_available, repo_root, test_options):
        """Test that filter_size_srna_sample is in the DAG for sized bigwigs."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_21NT_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "filter_size_srna_sample" in output, \
            "filter_size_srna_sample should be in the DAG for sized bigwigs"


# ===========================================================================
# TestReplicateMerging
# ===========================================================================

class TestReplicateMerging:
    """Test replicate merging rules for sRNA samples with multiple replicates."""

    def test_wt_two_reps_trigger_merging(self, snakemake_available, repo_root, test_options):
        """WT has 2 reps -- merged replicate rules should appear."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_MERGED_21NT_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "merging_srna_replicates" in output, \
            "sRNA with 2 reps should trigger merging_srna_replicates rule"

    def test_mutant_two_reps_trigger_merging(self, snakemake_available, repo_root, test_options):
        """mutant has 2 reps -- merged replicate rules should appear."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, MUTANT_MERGED_21NT_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "merging_srna_replicates" in output, \
            "sRNA with 2 reps should trigger merging_srna_replicates rule"


# ===========================================================================
# TestWildcardResolution
# ===========================================================================

class TestWildcardResolution:
    """Test wildcard resolution for sRNA samples."""

    def test_wt_rep1_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for WT rep1 at various sizes."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/sRNA/tracks/WT_sRNA_rep1__21nt__plus.bw",
            f"{_OUTPUT_DIR}/sRNA/tracks/WT_sRNA_rep1__22nt__minus.bw",
            f"{_OUTPUT_DIR}/sRNA/tracks/WT_sRNA_rep1__24nt__plus.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_mutant_rep1_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for mutant rep1."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/sRNA/tracks/mutant_sRNA_rep1__21nt__plus.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode == 0, f"Wildcard resolution failed for mutant: {result.stderr}"

    def test_mutant2_rep1_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for mutant2 rep1."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/sRNA/tracks/mutant2_sRNA_rep1__21nt__plus.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode == 0, f"Wildcard resolution failed for mutant2: {result.stderr}"

    def test_size_stats_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that size stats target resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_SIZE_STATS)
        assert result.returncode == 0, f"Wildcard resolution failed for size stats: {result.stderr}"

    def test_cluster_bed_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that cluster BED target resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_CLUSTERS)
        assert result.returncode == 0, f"Wildcard resolution failed for cluster BED: {result.stderr}"


# ===========================================================================
# TestDAGStructure
# ===========================================================================

class TestDAGStructure:
    """Test DAG generation and structure for sRNA workflow."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WT_REP1_21NT_PLUS)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_srna_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG contains sRNA-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WT_REP1_21NT_PLUS)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        for rule in ("shortstack_map", "filter_size_srna_sample", "make_srna_stranded_bigwigs"):
            assert rule in dag_output, f"Rule '{rule}' not found in sRNA DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_options):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WT_REP1_21NT_PLUS)
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
        target = f"{_OUTPUT_DIR}/sRNA/tracks/nonexistent_sample__21nt__plus.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for non-existent sample"

    def test_invalid_strand(self, snakemake_available, repo_root, test_options):
        """Test that an invalid strand label fails."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/sRNA/tracks/WT_sRNA_rep1__21nt__invalid_strand.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for invalid strand"


class TestShortStackOutputNamesMatchAligner:
    """ShortStack names its outputs from the genome-free input FASTQ.

    The rule declares them genome-qualified, so the shell must rename across
    that gap; a dry-run cannot run ShortStack, so check the commands (#71).
    """

    MV_RE = re.compile(r'mv "([^"]+)" "([^"]+)"')

    @pytest.fixture(scope="class")
    def shell_text(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(
            repo_root, test_options, WT_REP1_CLUSTERS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        return result.stdout + result.stderr

    def _condensed_moves(self, text):
        return [(a, b) for a, b in self.MV_RE.findall(text)
                if a.endswith("_condensed.bam") or a.endswith("_condensed.bam.csi")]

    def test_condensed_bam_is_renamed(self, shell_text):
        moves = self._condensed_moves(shell_text)
        assert moves, "no rename of the ShortStack condensed BAM found"
        for src, dst in moves:
            assert src != dst, f"pointless rename: {src}"

    def test_rename_source_carries_no_genome(self, shell_text):
        moves = self._condensed_moves(shell_text)
        assert moves, "no rename of the ShortStack condensed BAM found"
        for src, _ in moves:
            assert "__test_genome" not in os.path.basename(src), (
                f"rename source {src} is genome-qualified; ShortStack will not "
                "write that name")

    def test_rename_target_is_genome_qualified(self, shell_text):
        moves = self._condensed_moves(shell_text)
        assert moves, "no rename of the ShortStack condensed BAM found"
        for _, dst in moves:
            assert "__test_genome" in os.path.basename(dst), (
                f"rename target {dst} is not genome-qualified")

    def test_rename_source_matches_the_input_fastq_stem(self, shell_text):
        moves = self._condensed_moves(shell_text)
        assert moves, "no rename of the ShortStack condensed BAM found"
        fastqs = re.findall(r"(\S+/sRNA/fastq/clean__\S+?\.fastq\.gz)", shell_text)
        assert fastqs, "no clean FASTQ found in the shell commands"
        stems = {os.path.basename(f)[: -len(".fastq.gz")] for f in fastqs}
        for src, _ in moves:
            base = os.path.basename(src)
            assert any(base.startswith(stem) for stem in stems), (
                f"{base} does not start with any input FASTQ stem {sorted(stems)}")
