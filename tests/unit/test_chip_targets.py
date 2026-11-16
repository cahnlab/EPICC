"""Regression tests for ChIP peak-target selection.

Covers the "No control found" DAG-build failure: control rows that no IP
references ("orphan" controls) — and IPs missing a Control — were enumerated as
peak-calling targets, then crashed in assign_chip_input because a control row
has no Control of its own.
"""

import os
import sys

import pandas as pd
import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))
# samplefile_validation imports `from scripts.sample_sheet import ...`, so the
# workflow/ dir must also be importable (mirrors test_samplefile_validation.py).
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))

from sample_sheet import (  # noqa: E402
    add_compat_columns,
    build_analysis_name,
    get_analysis_samples,
    get_control_sample_id,
    get_replicate_sample_ids,
    identify_control_samples,
    is_peak_call_target,
    peak_callable_rows,
    read_sample_sheet,
)
from scripts.samplefile_validation import check_table  # noqa: E402
_CHIP_SHEET = os.path.join(
    _REPO_ROOT, "tests", "integration", "data", "test_samples_ChIP.tsv"
)


def _row(**kw):
    # Post-add_compat_columns shape: Assay carries the *combined* token, which is
    # what is_peak_call_target actually receives at runtime. Sheet-level fixtures
    # (which go through check_table) use the separated Assay + Peak_type form.
    base = {"Sample_ID": "s", "Assay": "ChIP_broad", "IP_target": "", "Control": ""}
    base.update(kw)
    return pd.Series(base)


class TestIsPeakCallTarget:
    def test_ip_with_control_is_target(self):
        assert is_peak_call_target(
            _row(Assay="ChIP_broad", IP_target="H3K9me2", Control="IN1")
        )

    def test_orphan_control_is_not_target(self):
        # A control row: no Control of its own -> cannot be peak-called.
        assert not is_peak_call_target(
            _row(Assay="ChIP_broad", IP_target="Input", Control="")
        )

    def test_ip_without_control_is_not_target(self):
        # An IP that forgot its Control would crash the same way.
        assert not is_peak_call_target(
            _row(Assay="ChIP_narrow", IP_target="H3K4me3", Control="")
        )

    def test_nan_control_is_not_target(self):
        assert not is_peak_call_target(
            _row(Assay="CUT_RUN_broad", IP_target="IgG", Control=float("nan"))
        )

    @pytest.mark.parametrize("assay", ["ATAC", "RNAseq", "sRNA", "WGBS", "dmC"])
    def test_non_pulldown_assays_always_targets(self, assay):
        # The Control column does not apply to these; they must not be filtered
        # out (ATAC in particular calls peaks with no control).
        assert is_peak_call_target(_row(Assay=assay, IP_target="", Control=""))


class TestChipSheetPeakTargets:
    """End-to-end on the real ChIP fixture, which contains two orphan controls."""

    @pytest.fixture
    def df(self):
        return add_compat_columns(read_sample_sheet(_CHIP_SHEET))

    def test_fixture_has_orphan_controls(self, df):
        referenced = identify_control_samples(df)
        control_like = df[
            df["IP_target"].str.strip().str.lower().isin(["input", "wce", "igg"])
        ]["Sample_ID"]
        orphans = sorted(s for s in control_like if s not in referenced)
        assert orphans == ["mutant_Input_rep1", "mutant_WCE_rep2"]

    def test_orphan_controls_excluded_from_peak_targets(self, df):
        controls = identify_control_samples(df)
        mask = (~df["Sample_ID"].isin(controls)) & df.apply(
            is_peak_call_target, axis=1
        )
        targets = set(df[mask]["Sample_ID"])
        assert "mutant_Input_rep1" not in targets
        assert "mutant_WCE_rep2" not in targets

    def test_real_ips_still_peak_targets(self, df):
        controls = identify_control_samples(df)
        mask = (~df["Sample_ID"].isin(controls)) & df.apply(
            is_peak_call_target, axis=1
        )
        targets = set(df[mask]["Sample_ID"])
        for sid in (
            "WT_H3K9me2_rep1",
            "WT_H3K9me2_rep2",
            "WT_H3K4me3_rep1",
            "mutant_H3K9me2_rep1",
            "mutant_H3K4me3_rep1",
        ):
            assert sid in targets, f"{sid} must remain a peak target"

    def test_orphan_controls_excluded_from_analysis_samples(self, df):
        # They also leaked in here, producing selected_peaks/IDR targets.
        names = set(get_analysis_samples(df)["sample_name"])
        assert "mutant_Input_rep1" not in names
        assert "mutant_WCE_rep2" not in names
        assert "WT_H3K9me2_rep1" in names


class TestNonPulldownAnalysisSamplesUnaffected:
    """get_analysis_samples is shared by ATAC/RNA/sRNA/mC — those envs have no
    Control column semantics and must be untouched by the new filter."""

    def test_rna_and_mc_samples_retained(self):
        # Distinct Levels per row so replicate dedup doesn't collapse them.
        df = pd.DataFrame(
            {
                "Sample_ID": ["r1", "r2", "m1", "a1"],
                "Assay": ["RNAseq", "RNAseq", "WGBS", "ATAC"],
                "Genome": ["G"] * 4,
                "Levels": [
                    "genotype:WT",
                    "genotype:mut",
                    "genotype:WT2",
                    "genotype:WT3",
                ],
                "Replicate_ID": ["rep1"] * 4,
                "Read_files": ["SRR1", "SRR2", "SRR3", "SRR4"],
                "Read_layout": ["SE"] * 4,
                "IP_target": [""] * 4,
                "Control": [""] * 4,
            }
        )
        names = set(get_analysis_samples(add_compat_columns(df))["Sample_ID"])
        assert names == {"r1", "r2", "m1", "a1"}

    def test_new_filter_is_noop_for_non_pulldown_sheets(self):
        """The added is_peak_call_target filter must not change which analysis
        samples a control-free (non-pulldown) sheet yields."""
        df = add_compat_columns(
            pd.DataFrame(
                {
                    "Sample_ID": ["r1", "r2", "m1"],
                    "Assay": ["RNAseq", "RNAseq", "EMseq"],
                    "Genome": ["G"] * 3,
                    "Levels": ["genotype:WT", "genotype:mut", "genotype:WT"],
                    "Replicate_ID": ["rep1"] * 3,
                    "Read_files": ["SRR1", "SRR2", "SRR3"],
                    "Read_layout": ["SE"] * 3,
                    "IP_target": [""] * 3,
                    "Control": [""] * 3,
                }
            )
        )
        # Every row passes the predicate, so nothing is filtered.
        assert df.apply(is_peak_call_target, axis=1).all()
        assert len(get_analysis_samples(df)) == 3


class TestNoControlWarning:
    """check_table warns (not errors) when a pulldown IP has no Control, since
    such a sample is silently dropped from the peak-target set."""

    @staticmethod
    def _sheet(*specs):
        return pd.DataFrame(
            [
                {
                    "Sample_ID": sid,
                    "Assay": assay,
                    "Genome": "G",
                    "Levels": "genotype:WT",
                    "Replicate_ID": "rep1",
                    "Read_files": f"SRR{n + 1}",
                    "Read_layout": "SE",
                    "IP_target": ip,
                    "Control": ctrl,
                    "Peak_type": "broad",
                }
                for n, (sid, assay, ip, ctrl) in enumerate(specs)
            ]
        )

    def _warnings(self, df, capsys):
        check_table(df, check_paths=False)  # must not raise
        return [l for l in capsys.readouterr().out.splitlines() if l.startswith("[!]")]

    def test_warns_for_ip_without_control(self, capsys):
        df = self._sheet(("IP1", "ChIP", "H3K9me2", ""))
        warns = self._warnings(df, capsys)
        assert len(warns) == 1
        assert "will not be peak-called" in warns[0]
        assert "IP1" in warns[0]

    @pytest.mark.parametrize("ip_target", ["Input", "WCE", "IgG", "input", "igg"])
    def test_no_warning_for_control_rows(self, ip_target, capsys):
        # Control rows legitimately have no Control of their own.
        df = self._sheet(("C1", "ChIP", ip_target, ""))
        assert self._warnings(df, capsys) == []

    def test_no_warning_for_proper_ip_control_pair(self, capsys):
        df = self._sheet(
            ("IP1", "ChIP", "H3K9me2", "IN1"),
            ("IN1", "ChIP", "Input", ""),
        )
        assert self._warnings(df, capsys) == []

    def test_no_warning_on_real_chip_fixture(self, capsys):
        # Every IP there has a Control; the orphans are Input/WCE rows.
        df = read_sample_sheet(_CHIP_SHEET)
        assert self._warnings(df, capsys) == []


class TestDualRoleControlAndTarget:
    """A sample may serve as another row's control AND be analysed itself.

    Real case: an H3 ChIP used both to see raw H3 distribution (its own peaks,
    tracks, merged analysis, IDR) and as the peak-calling control for H3K9me2.
    Qualifying is structural — the row must declare a ``Control`` of its own —
    so being referenced by someone else is not disqualifying.
    """

    HEADER = ("Assay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
              "Read_layout\tIP_target\tControl")
    ROWS = [
        # Input: no Control of its own -> not analysable
        "ChIP_broad\tG\tgenotype:WT\trep1\tSRR1\tSE\tInput\t",
        # H3: dual role -> control for H3K9me2, but has its own Control (Input)
        "ChIP_broad\tG\tgenotype:WT\trep1\tSRR2\tSE\tH3\tWT_Input_rep1",
        "ChIP_broad\tG\tgenotype:WT\trep2\tSRR3\tSE\tH3\tWT_Input_rep1",
        # the IP that uses H3 as its control
        "ChIP_broad\tG\tgenotype:WT\trep1\tSRR4\tSE\tH3K9me2\tWT_H3_rep1",
        # orphan control (regression guard for the #53 fix)
        "ChIP_broad\tG\tgenotype:WT\trep1\tSRR5\tSE\tWCE\t",
    ]

    @pytest.fixture
    def df(self, tmp_path):
        f = tmp_path / "dual.tsv"
        f.write_text("\n".join([self.HEADER] + self.ROWS) + "\n")
        return add_compat_columns(read_sample_sheet(f))

    def test_dual_role_sample_is_referenced_as_a_control(self, df):
        # Precondition: WT_H3_rep1 really is someone's control.
        assert "WT_H3_rep1" in identify_control_samples(df)

    def test_dual_role_reps_get_peak_targets(self, df):
        targets = set(peak_callable_rows(df)["Sample_ID"])
        assert {"WT_H3_rep1", "WT_H3_rep2"} <= targets, \
            "dual-role sample must get its own peak/FC targets"
        assert "WT_H3K9me2_rep1" in targets

    def test_unanalysable_rows_still_excluded(self, df):
        targets = set(peak_callable_rows(df)["Sample_ID"])
        # Input has no Control; the orphan WCE is referenced by nobody and has
        # no Control -- neither can be peak-called.
        assert "WT_Input_rep1" not in targets
        assert "WT_WCE_rep1" not in targets

    def test_dual_role_gets_analysis_level_treatment(self, df):
        analysis = get_analysis_samples(df)
        names = {build_analysis_name(r) for _, r in analysis.iterrows()}
        assert "ChIP_broad__WT__H3__G" in names, \
            "dual-role sample must get a merged analysis group (peaks/IDR/plots)"
        assert "ChIP_broad__WT__H3K9me2__G" in names

    def test_dual_role_analysis_group_has_both_replicates(self, df):
        # Merged analysis / IDR need both H3 reps grouped under the analysis name.
        # The returned names are genome-qualified (mapped_name): every caller
        # uses them to build post-alignment paths, where the reference is part
        # of the file's identity.
        assert sorted(get_replicate_sample_ids("ChIP_broad__WT__H3__G", df)) == \
            ["WT_H3_rep1__G", "WT_H3_rep2__G"]

    def test_dual_role_still_resolvable_as_a_control(self, df):
        # Looking the sample up by Sample_ID must still hit the control-merge
        # path, so H3K9me2's control resolution is unaffected. Sample_IDs cannot
        # contain '__' while analysis names always do, so the two lookup key
        # formats can never collide. The *key* stays bare; the returned names
        # are genome-qualified.
        assert get_replicate_sample_ids("WT_H3_rep1", df) == ["WT_H3_rep1__G"]

    def test_dual_role_control_resolution_stays_genome_matched(self, df):
        # The qualified name a caller got from get_replicate_sample_ids must
        # round-trip through control lookup with the genome preserved — a
        # control aligned to a different reference would be the wrong BAM.
        assert get_control_sample_id("WT_H3K9me2_rep1__G", df) == "WT_H3_rep1__G"
        assert get_control_sample_id("WT_H3_rep1__G", df) == "WT_Input_rep1__G"
