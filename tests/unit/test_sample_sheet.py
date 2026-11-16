"""Unit tests for workflow/scripts/sample_sheet.py"""

import pytest
import pandas as pd
from collections import OrderedDict
import sys
import os

# Add the workflow/scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "scripts"))
# samplefile_validation imports `from scripts.sample_sheet import ...`, which
# only resolves with the parent of the workflow/ tree on sys.path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow"))

from scripts.samplefile_validation import check_table  # noqa: E402

from sample_sheet import (
    VALID_ASSAYS,
    ASSAY_TO_ENV,
    ASSAY_TO_PEAKTYPE,
    NEW_COLNAMES,
    parse_levels,
    levels_to_label,
    levels_to_factors,
    parse_read_files,
    get_seq_id_and_path,
    build_analysis_key,
    build_analysis_name,
    build_control_merge_key,
    build_analysis_to_replicates,
    get_replicate_sample_ids,
    identify_control_samples,
    get_control_sample_id,
    get_analysis_samples,
    get_sample_field,
    get_env,
    get_peaktype,
    add_compat_columns,
    read_sample_sheet,
    s3_uri_to_https,
    combine_assay_peaktype,
    PEAK_TYPE_ASSAYS,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def pombe_df():
    """A minimal pombe-like sample sheet in new format."""
    data = {
        "Sample_ID": [
            "WT_H3K9me2_rep1", "WT_H3K9me2_rep2",
            "WT_WCE_rep1",
            "dcr1_H3K9me2_rep1",
            "WT_RNA_rep1", "WT_RNA_rep2",
            "dcr1_RNA_rep1",
            "WT_sRNA_rep1",
        ],
        "Assay": [
            "ChIP_broad", "ChIP_broad",
            "ChIP_broad",
            "ChIP_broad",
            "RNAseq", "RNAseq",
            "RNAseq",
            "sRNA",
        ],
        "Genome": ["Spombe"] * 8,
        "Levels": [
            "genotype:WT", "genotype:WT",
            "genotype:WT",
            "genotype:dcr1",
            "genotype:WT", "genotype:WT",
            "genotype:dcr1",
            "genotype:WT",
        ],
        "Replicate_ID": [
            "rep1", "rep2", "rep1", "rep1",
            "rep1", "rep2", "rep1", "rep1",
        ],
        "Read_files": [
            "SRR20678305", "SRR20678333",
            "SRR5445712",
            "SRR20678308",
            "SRR30889044", "SRR30889043",
            "SRR30889039",
            "SRR20678362",
        ],
        "Read_layout": [
            "PE", "PE", "PE", "PE",
            "SE", "SE", "SE", "SE",
        ],
        "IP_target": [
            "H3K9me2", "H3K9me2", "WCE",
            "H3K9me2",
            "", "", "", "",
        ],
        "Control": [
            "WT_WCE_rep1", "WT_WCE_rep1", "",
            "WT_WCE_rep1",
            "", "", "", "",
        ],
    }
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# parse_levels
# ---------------------------------------------------------------------------

class TestParseLevels:
    def test_basic(self):
        result = parse_levels("genotype:WT,tissue:root")
        assert result == OrderedDict([("genotype", "WT"), ("tissue", "root")])

    def test_single_factor(self):
        result = parse_levels("genotype:WT")
        assert result == OrderedDict([("genotype", "WT")])

    def test_three_factors(self):
        result = parse_levels("genotype:WT,tissue:root,time:T0")
        assert len(result) == 3
        assert result["time"] == "T0"

    def test_whitespace_tolerant(self):
        result = parse_levels("genotype: WT , tissue: root")
        assert result == OrderedDict([("genotype", "WT"), ("tissue", "root")])

    def test_invalid_no_colon(self):
        with pytest.raises(ValueError, match="expected 'factor:level'"):
            parse_levels("genotype_WT")


class TestLevelsToLabel:
    def test_basic(self):
        assert levels_to_label("genotype:WT,tissue:root") == "WT_root"

    def test_single(self):
        assert levels_to_label("genotype:WT") == "WT"

    def test_three(self):
        assert levels_to_label("genotype:WT,tissue:root,time:T0") == "WT_root_T0"


class TestLevelsToFactors:
    def test_basic(self):
        assert levels_to_factors("genotype:WT,tissue:root") == ["genotype", "tissue"]


# ---------------------------------------------------------------------------
# parse_read_files
# ---------------------------------------------------------------------------

class TestParseReadFiles:
    def test_single_sra(self):
        parts, is_sra = parse_read_files("SRR20678305", "PE")
        assert parts == ["SRR20678305"]
        assert is_sra is True

    def test_merged_sra(self):
        parts, is_sra = parse_read_files("SRR111+SRR222+SRR333", "SE")
        assert parts == ["SRR111", "SRR222", "SRR333"]
        assert is_sra is True

    def test_ddbj_accession(self):
        parts, is_sra = parse_read_files("DRR400324", "PE")
        assert parts == ["DRR400324"]
        assert is_sra is True

    def test_ena_accession(self):
        parts, is_sra = parse_read_files("ERR123456", "SE")
        assert parts == ["ERR123456"]
        assert is_sra is True

    def test_local_se(self):
        parts, is_sra = parse_read_files("/path/to/reads.fq.gz", "SE")
        assert parts == ["/path/to/reads.fq.gz"]
        assert is_sra is False

    def test_local_pe(self):
        parts, is_sra = parse_read_files("/path/r1.fq.gz,/path/r2.fq.gz", "PE")
        assert parts == ["/path/r1.fq.gz,/path/r2.fq.gz"]
        assert is_sra is False


class TestGetSeqIdAndPath:
    def test_sra(self):
        seq_id, fq_path = get_seq_id_and_path("SRR20678305", "PE")
        assert seq_id == "SRR20678305"
        assert fq_path == "SRA"

    def test_sra_merged(self):
        seq_id, fq_path = get_seq_id_and_path("SRR111+SRR222", "SE")
        assert seq_id == "SRR111,SRR222"
        assert fq_path == "SRA"

    def test_ddbj(self):
        seq_id, fq_path = get_seq_id_and_path("DRR400324", "PE")
        assert seq_id == "DRR400324"
        assert fq_path == "SRA"

    def test_bam(self):
        seq_id, fq_path = get_seq_id_and_path("/path/to/sample.bam", "SE")
        assert seq_id == "sample"
        assert fq_path == "/path/to/sample.bam"

    def test_bam_url(self):
        url = "https://lemna.org/data/shimada2024/SRR28453434.bam"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "SRR28453434"
        assert fq_path == url

    def test_bedmethyl_url(self):
        url = "https://example.com/data/sample.bed.gz"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "sample"
        assert fq_path == url

    def test_bedmethyl_url_with_query_params(self):
        url = "https://example.com/data/sample.bed.gz?token=abc123"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "sample"
        assert fq_path == url

    def test_fastq_url_se(self):
        url = "https://example.com/reads.fastq.gz"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "URL"
        assert fq_path == url

    def test_fastq_url_pe(self):
        r1 = "https://example.com/r1.fq.gz"
        r2 = "https://example.com/r2.fq.gz"
        seq_id, fq_path = get_seq_id_and_path(f"{r1},{r2}", "PE")
        assert seq_id == "URL"
        assert fq_path == f"{r1},{r2}"


# ---------------------------------------------------------------------------
# Analysis key / name
# ---------------------------------------------------------------------------

class TestBuildAnalysisKey:
    def test_chip(self, pombe_df):
        row = pombe_df.iloc[0]
        key = build_analysis_key(row)
        assert key == ("ChIP_broad", "WT", "H3K9me2", "Spombe")

    def test_rnaseq(self, pombe_df):
        row = pombe_df.iloc[5]
        key = build_analysis_key(row)
        assert key == ("RNAseq", "WT", "", "Spombe")


class TestBuildAnalysisName:
    def test_chip(self, pombe_df):
        row = pombe_df.iloc[0]
        name = build_analysis_name(row)
        assert name == "ChIP_broad__WT__H3K9me2__Spombe"

    def test_rnaseq(self, pombe_df):
        row = pombe_df.iloc[5]
        name = build_analysis_name(row)
        assert name == "RNAseq__WT__Spombe"


# ---------------------------------------------------------------------------
# Analysis-to-replicates
# ---------------------------------------------------------------------------

class TestBuildAnalysisToReplicates:
    def test_basic(self, pombe_df):
        a2r = build_analysis_to_replicates(pombe_df)
        # WT H3K9me2 ChIP has 2 reps
        key = ("ChIP_broad", "WT", "H3K9me2", "Spombe")
        assert key in a2r
        assert a2r[key] == ["rep1", "rep2"]

    def test_controls_excluded(self, pombe_df):
        a2r = build_analysis_to_replicates(pombe_df)
        # Input samples should not appear as analysis keys
        for key in a2r:
            assert key[2] != "WCE"  # IP_target should never be a control type

    def test_single_rep(self, pombe_df):
        a2r = build_analysis_to_replicates(pombe_df)
        key = ("ChIP_broad", "dcr1", "H3K9me2", "Spombe")
        assert key in a2r
        assert a2r[key] == ["rep1"]


# ---------------------------------------------------------------------------
# Control handling
# ---------------------------------------------------------------------------

class TestIdentifyControlSamples:
    def test_basic(self, pombe_df):
        controls = identify_control_samples(pombe_df)
        assert controls == {"WT_WCE_rep1"}


class TestGetControlSampleId:
    def test_has_control(self, pombe_df):
        ctrl = get_control_sample_id("WT_H3K9me2_rep1", pombe_df)
        assert ctrl == "WT_WCE_rep1"

    def test_no_control(self, pombe_df):
        ctrl = get_control_sample_id("WT_RNA_rep1", pombe_df)
        assert ctrl is None

    def test_is_control(self, pombe_df):
        ctrl = get_control_sample_id("WT_WCE_rep1", pombe_df)
        assert ctrl is None


class TestControlMergeAcrossAssay:
    """Control rows that represent the same biological material (Input/WCE/IgG)
    should merge into one replicate group even if individual rep rows are
    labeled with different Assay values (e.g. one ChIP_broad, one ChIP_narrow).
    See build_control_merge_key — Assay is intentionally not part of the key.
    """

    @pytest.fixture
    def mixed_assay_df(self):
        # rep1 Input labeled ChIP_broad, rep2 Input labeled ChIP_narrow —
        # biologically the same Input pool. Two IPs (one of each Assay) each
        # point at one of the Input reps via the Control column.
        return pd.DataFrame({
            "Sample_ID": [
                "WT_leaf_H3K9me2_rep1", "WT_leaf_H3K4me3_rep1",
                "WT_leaf_Input_rep1", "WT_leaf_Input_rep2",
            ],
            "Assay": [
                "ChIP_broad", "ChIP_narrow",
                "ChIP_broad", "ChIP_narrow",
            ],
            "Genome": ["ColCEN"] * 4,
            "Levels": ["genotype:WT,tissue:leaf"] * 4,
            "Replicate_ID": ["rep1", "rep1", "rep1", "rep2"],
            "Read_files": ["SRR1", "SRR2", "SRR3", "SRR4"],
            "Read_layout": ["PE"] * 4,
            "IP_target": ["H3K9me2", "H3K4me3", "Input", "Input"],
            "Control": [
                "WT_leaf_Input_rep1", "WT_leaf_Input_rep2",
                "", "",
            ],
        })

    def test_merge_key_drops_assay(self, mixed_assay_df):
        rep1_key = build_control_merge_key(mixed_assay_df.iloc[2])
        rep2_key = build_control_merge_key(mixed_assay_df.iloc[3])
        assert rep1_key == rep2_key == ("WT_leaf", "Input", "ColCEN")

    def test_merge_key_differs_from_full_analysis_key(self, mixed_assay_df):
        # build_analysis_key still distinguishes Assay — only the
        # control-merge variant strips it.
        rep1_full = build_analysis_key(mixed_assay_df.iloc[2])
        rep2_full = build_analysis_key(mixed_assay_df.iloc[3])
        assert rep1_full != rep2_full

    def test_get_replicate_sample_ids_merges_across_assay(self, mixed_assay_df):
        # Asking for replicates of EITHER Input rep returns BOTH — they
        # share the biological-control merge group.
        for ctrl in ("WT_leaf_Input_rep1", "WT_leaf_Input_rep2"):
            reps = get_replicate_sample_ids(ctrl, mixed_assay_df)
            assert set(reps) == {"WT_leaf_Input_rep1", "WT_leaf_Input_rep2"}, (
                f"Expected both Input reps in the merge group for {ctrl}, got {reps}"
            )

    def test_different_ip_target_does_not_merge(self, mixed_assay_df):
        # A separate IgG control with the same Levels/Genome but different
        # IP_target must NOT merge with the Input group.
        igg_row = pd.DataFrame({
            "Sample_ID": ["WT_leaf_IgG_rep1"], "Assay": ["ChIP_broad"],
            "Genome": ["ColCEN"], "Levels": ["genotype:WT,tissue:leaf"],
            "Replicate_ID": ["rep1"], "Read_files": ["SRR5"],
            "Read_layout": ["PE"], "IP_target": ["IgG"],
            "Control": [""],
        })
        df_extended = pd.concat([mixed_assay_df, igg_row], ignore_index=True)
        # Add an IP that references the IgG so it shows up in identify_control_samples
        df_extended = pd.concat([df_extended, pd.DataFrame({
            "Sample_ID": ["WT_leaf_TF_rep1"], "Assay": ["ChIP_narrow"],
            "Genome": ["ColCEN"], "Levels": ["genotype:WT,tissue:leaf"],
            "Replicate_ID": ["rep1"], "Read_files": ["SRR6"],
            "Read_layout": ["PE"], "IP_target": ["TF1"],
            "Control": ["WT_leaf_IgG_rep1"],
        })], ignore_index=True)

        input_reps = get_replicate_sample_ids("WT_leaf_Input_rep1", df_extended)
        igg_reps = get_replicate_sample_ids("WT_leaf_IgG_rep1", df_extended)
        assert "WT_leaf_IgG_rep1" not in input_reps
        assert "WT_leaf_Input_rep1" not in igg_reps


# ---------------------------------------------------------------------------
# Analysis samples (non-control)
# ---------------------------------------------------------------------------

class TestGetAnalysisSamples:
    def test_controls_excluded(self, pombe_df):
        analysis = get_analysis_samples(pombe_df)
        # Controls should not appear
        assert "WT_WCE_rep1" not in analysis["Sample_ID"].values
        # dcr1_WCE_rep1 removed — no dcr1 WCE in public data; dcr1 samples use WT_WCE_rep1

    def test_deduplicated(self, pombe_df):
        analysis = get_analysis_samples(pombe_df)
        # WT H3K9me2 has 2 reps but should appear once in analysis
        chip_wt = analysis[
            (analysis["Assay"] == "ChIP_broad") &
            (analysis["IP_target"] == "H3K9me2") &
            (analysis["Levels"].str.contains("WT"))
        ]
        assert len(chip_wt) == 1


# ---------------------------------------------------------------------------
# Field lookup
# ---------------------------------------------------------------------------

class TestGetSampleField:
    def test_basic(self, pombe_df):
        val = get_sample_field("WT_H3K9me2_rep1", pombe_df, "Assay")
        assert val == "ChIP_broad"

    def test_missing(self, pombe_df):
        val = get_sample_field("nonexistent", pombe_df, "Assay")
        assert val is None


# ---------------------------------------------------------------------------
# Environment and peak type
# ---------------------------------------------------------------------------

class TestGetEnv:
    def test_chip_broad(self):
        assert get_env("ChIP_broad") == "ChIP"

    def test_rnaseq(self):
        assert get_env("RNAseq") == "RNA"

    def test_unknown(self):
        with pytest.raises(ValueError, match="Unknown assay"):
            get_env("unknown_assay")


class TestGetPeaktype:
    def test_chip_broad(self):
        assert get_peaktype("ChIP_broad") == "broad"

    def test_chip_narrow(self):
        assert get_peaktype("ChIP_narrow") == "narrow"

    def test_atac(self):
        assert get_peaktype("ATAC") == "narrow"

    def test_cut_run_broad(self):
        assert get_peaktype("CUT_RUN_broad") == "broad"

    def test_cut_run_narrow(self):
        assert get_peaktype("CUT_RUN_narrow") == "narrow"

    def test_cut_tag_broad(self):
        assert get_peaktype("CUT_TAG_broad") == "broad"

    def test_cut_tag_narrow(self):
        assert get_peaktype("CUT_TAG_narrow") == "narrow"

    def test_config_override(self):
        assert get_peaktype("ChIP_broad", {"ChIP_broad": "narrow"}) == "narrow"

    def test_no_peaktype(self):
        with pytest.raises(ValueError, match="No peak type"):
            get_peaktype("RNAseq")


class TestCombinedAssayTokensRejectedOnInput:
    """The combined token is internal only — a sheet must use Assay + Peak_type.

    It is still the key that output paths and env/peaktype lookups are built
    from, so these tests pin the boundary: rejected as input, required
    internally.
    """

    HEADER = ("Assay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
              "Read_layout\tIP_target\tControl\tPeak_type")

    def _sheet(self, tmp_path, assay, peak_type):
        f = tmp_path / "s.tsv"
        f.write_text("\n".join([
            self.HEADER,
            f"{assay}\tG\tgenotype:WT\trep1\tSRR1\tSE\tInput\t\t{peak_type}",
            f"{assay}\tG\tgenotype:WT\trep1\tSRR2\tSE\tH3K9me2\tWT_ChIP_Input_rep1\t{peak_type}",
        ]) + "\n")
        return read_sample_sheet(f)

    @pytest.mark.parametrize("assay", [
        "ChIP_broad", "ChIP_narrow", "CUT_RUN_broad", "CUT_TAG_narrow",
    ])
    def test_combined_assay_is_a_validation_error(self, tmp_path, assay):
        df = self._sheet(tmp_path, assay, "")
        with pytest.raises(ValueError, match="not in"):
            check_table(df, check_paths=False)

    def test_separated_form_passes(self, tmp_path):
        df = self._sheet(tmp_path, "ChIP", "broad")
        check_table(df, check_paths=False)  # no raise

    def test_combine_still_produces_the_internal_token(self, tmp_path):
        """Rejecting it as input must not stop add_compat_columns making it."""
        df = add_compat_columns(self._sheet(tmp_path, "ChIP", "broad"))
        assert set(df["Assay"]) == {"ChIP_broad"}
        assert set(df["env"]) == {"ChIP"}


class TestCUTAssayVocabulary:
    """Sample_sheet should treat CUT&RUN/CUT&Tag as IP_PEAK_ASSAYS routed to the ChIP env."""

    def test_both_in_valid_assays(self):
        for a in ("CUT_RUN", "CUT_TAG"):
            assert a in VALID_ASSAYS

    def test_combined_tokens_are_not_accepted_input(self):
        """The combined form is internal only; a sheet must use Assay+Peak_type."""
        for a in ("ChIP_broad", "ChIP_narrow",
                  "CUT_RUN_broad", "CUT_RUN_narrow",
                  "CUT_TAG_broad", "CUT_TAG_narrow"):
            assert a not in VALID_ASSAYS

    def test_all_four_route_to_chip_env(self):
        for a in ("CUT_RUN_broad", "CUT_RUN_narrow",
                  "CUT_TAG_broad", "CUT_TAG_narrow"):
            assert ASSAY_TO_ENV[a] == "ChIP"

    def test_peaktype_suffix_matches_designation(self):
        for a in ("CUT_RUN_broad", "CUT_TAG_broad"):
            assert ASSAY_TO_PEAKTYPE[a] == "broad"
        for a in ("CUT_RUN_narrow", "CUT_TAG_narrow"):
            assert ASSAY_TO_PEAKTYPE[a] == "narrow"

    def test_in_ip_peak_assays_set(self):
        from workflow.scripts.sample_sheet import IP_PEAK_ASSAYS
        for a in ("ChIP_broad", "ChIP_narrow",
                  "CUT_RUN_broad", "CUT_RUN_narrow",
                  "CUT_TAG_broad", "CUT_TAG_narrow"):
            assert a in IP_PEAK_ASSAYS
        # Non-IP assays must not leak into the set
        for a in ("ATAC", "RNAseq", "RAMPAGE", "sRNA", "WGBS", "dmC"):
            assert a not in IP_PEAK_ASSAYS


# ---------------------------------------------------------------------------
# Compat columns
# ---------------------------------------------------------------------------

class TestAddCompatColumns:
    def test_basic_columns_exist(self, pombe_df):
        result = add_compat_columns(pombe_df)
        for col in ["data_type", "ref_genome", "replicate", "paired",
                     "line", "tissue", "sample_type", "extra_info",
                     "levels_label", "env", "sample_name",
                     "seq_id", "fastq_path"]:
            assert col in result.columns, f"Missing column: {col}"

    def test_sample_name_is_sample_id(self, pombe_df):
        result = add_compat_columns(pombe_df)
        assert (result["sample_name"] == result["Sample_ID"]).all()

    def test_env_mapping(self, pombe_df):
        result = add_compat_columns(pombe_df)
        chip_rows = result[result["Assay"] == "ChIP_broad"]
        assert (chip_rows["env"] == "ChIP").all()

    def test_line_from_levels(self, pombe_df):
        result = add_compat_columns(pombe_df)
        wt_row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert wt_row["line"] == "WT"
        assert wt_row["tissue"] == ""

    def test_sample_type_chip(self, pombe_df):
        result = add_compat_columns(pombe_df)
        wt_chip = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert wt_chip["sample_type"] == "H3K9me2"

    def test_sample_type_rna(self, pombe_df):
        result = add_compat_columns(pombe_df)
        rna_row = result[result["Sample_ID"] == "WT_RNA_rep1"].iloc[0]
        assert rna_row["sample_type"] == "RNAseq"

    def test_extra_info_chip(self, pombe_df):
        result = add_compat_columns(pombe_df)
        chip_row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert chip_row["extra_info"] == "H3K9me2"

    def test_extra_info_rna(self, pombe_df):
        result = add_compat_columns(pombe_df)
        rna_row = result[result["Sample_ID"] == "WT_RNA_rep1"].iloc[0]
        assert rna_row["extra_info"] == "N/A"

    def test_sra_seq_id(self, pombe_df):
        result = add_compat_columns(pombe_df)
        row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert row["seq_id"] == "SRR20678305"
        assert row["fastq_path"] == "SRA"

    def test_levels_label(self, pombe_df):
        result = add_compat_columns(pombe_df)
        row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert row["levels_label"] == "WT"


class TestReadSampleSheet:
    HEADER = "Assay\tGenome\tLevels\tReplicate_ID\tRead_files\tRead_layout\tIP_target\tControl"
    ROW_A = "ChIP_broad\tColCEN\tgenotype:WT\trep1\tSRR000001\tSE\tH3K9me2\tWT_ChIP_Input_rep1"
    ROW_B = "ChIP_broad\tColCEN\tgenotype:WT\trep2\tSRR000002\tSE\tH3K9me2\tWT_ChIP_Input_rep1"
    ROW_C = "ChIP_broad\tColCEN\tgenotype:WT\trep1\tSRR000003\tSE\tInput\t"

    def test_full_line_comments_skipped(self, tmp_path):
        f = tmp_path / "samples.tsv"
        f.write_text(
            "\n".join([
                "# this is a header comment",
                self.HEADER,
                "# parked sample, do not run",
                self.ROW_A,
                self.ROW_C,
            ]) + "\n"
        )
        df = read_sample_sheet(f)
        assert sorted(df["Sample_ID"]) == ["WT_ChIP_H3K9me2_rep1", "WT_ChIP_Input_rep1"]

    def test_no_comments_unaffected(self, tmp_path):
        f = tmp_path / "samples.tsv"
        f.write_text("\n".join([self.HEADER, self.ROW_A, self.ROW_B, self.ROW_C]) + "\n")
        df = read_sample_sheet(f)
        assert sorted(df["Sample_ID"]) == [
            "WT_ChIP_H3K9me2_rep1", "WT_ChIP_H3K9me2_rep2", "WT_ChIP_Input_rep1"]


class TestPeakTypeSeparation:
    """Separated Assay + Peak_type schema and its back-compat with the legacy
    combined form (Action item: separate peak-calling type from assay type)."""

    HEADER = ("Assay\tPeak_type\tGenome\tLevels\tReplicate_ID\t"
              "Read_files\tRead_layout\tIP_target\tControl")

    def _read(self, tmp_path, rows):
        f = tmp_path / "s.tsv"
        f.write_text("\n".join([self.HEADER] + rows) + "\n")
        return read_sample_sheet(f)

    def test_combine_helper(self):
        assert combine_assay_peaktype("ChIP", "broad") == "ChIP_broad"
        assert combine_assay_peaktype("CUT_TAG", "narrow") == "CUT_TAG_narrow"
        # idempotent on legacy combined tokens
        assert combine_assay_peaktype("ChIP_broad", "") == "ChIP_broad"
        # non-peak / ATAC / missing peak type: unchanged
        assert combine_assay_peaktype("ATAC", "") == "ATAC"
        assert combine_assay_peaktype("RNAseq", "") == "RNAseq"
        assert combine_assay_peaktype("ChIP", "") == "ChIP"

    def test_peak_type_assays_collapsed(self):
        assert PEAK_TYPE_ASSAYS == {"ChIP", "CUT_RUN", "CUT_TAG"}

    def test_separated_form_combines_and_keeps_peak_in_name(self, tmp_path):
        df = self._read(tmp_path, [
            "ChIP\tbroad\tG\tgenotype:WT\trep1\tSRR1\tSE\tH3K9me2\tWT_ChIP_Input_rep1",
            "ChIP\tbroad\tG\tgenotype:WT\trep1\tSRR2\tSE\tInput\t",
        ])
        c = add_compat_columns(df)
        assert list(c["Assay"]) == ["ChIP_broad", "ChIP_broad"]
        assert set(c["env"]) == {"ChIP"}
        ip = c[c["Sample_ID"] == "WT_ChIP_H3K9me2_rep1"].iloc[0]
        assert build_analysis_name(ip) == "ChIP_broad__WT__H3K9me2__G"

    def test_legacy_combined_still_works(self, tmp_path):
        # Legacy sheet: peak type baked into Assay, no Peak_type column.
        f = tmp_path / "legacy.tsv"
        f.write_text(
            "Assay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
            "Read_layout\tIP_target\tControl\n"
            "ChIP_narrow\tG\tgenotype:WT\trep1\tSRR1\tSE\tH3K4me3\tWT_ChIP_Input_rep1\n"
            "ChIP_narrow\tG\tgenotype:WT\trep1\tSRR2\tSE\tInput\t\n"
        )
        c = add_compat_columns(read_sample_sheet(f))
        assert list(c["Assay"]) == ["ChIP_narrow", "ChIP_narrow"]
        assert get_peaktype("ChIP_narrow") == "narrow"

    def test_new_form_peaktype_via_combined(self, tmp_path):
        df = self._read(tmp_path, [
            "CUT_TAG\tnarrow\tG\tgenotype:WT\trep1\tSRR1\tSE\tCTCF\t",
        ])
        c = add_compat_columns(df)
        assert get_peaktype(c.iloc[0]["Assay"]) == "narrow"
