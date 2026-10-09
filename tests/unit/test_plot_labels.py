"""Tests for plot_levels_labels / check_unique_labels in sample_sheet.py.

Combined-analysis labels are levels_label plus the mark, so two assays
profiling the same mark from the same Levels would share one label.
"""

import os
import sys

import pandas as pd
import pytest

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(REPO, "workflow", "scripts"))

from sample_sheet import (  # noqa: E402
    add_compat_columns, check_unique_labels, get_analysis_samples,
    plot_levels_labels, read_sample_sheet,
)

DATA = os.path.join(REPO, "tests", "integration", "data")


def _analysis(sheet):
    samples = add_compat_columns(read_sample_sheet(os.path.join(DATA, sheet)))
    return samples, add_compat_columns(get_analysis_samples(samples))


def _rows(*rows):
    return pd.DataFrame(rows, columns=["Assay", "env", "sample_type", "levels_label"])


class TestPlotLevelsLabels:
    def test_no_collision_unchanged(self):
        df = _rows(("ChIP_broad", "ChIP", "H3K9me2", "WT"),
                   ("ChIP_broad", "ChIP", "H3K9me2", "dcr1"),
                   ("ChIP_narrow", "ChIP", "H3K4me3", "WT"))
        assert list(plot_levels_labels(df)) == ["WT", "dcr1", "WT"]

    def test_pulldown_assays_get_assay_without_peak_type(self):
        df = _rows(("CUT_RUN_broad", "ChIP", "H3K27me3", "WT"),
                   ("CUT_TAG_broad", "ChIP", "H3K27me3", "WT"),
                   ("CUT_RUN_broad", "ChIP", "H3K27me3", "mut"))
        assert list(plot_levels_labels(df)) == ["WT_CUT_RUN", "WT_CUT_TAG", "mut"]

    def test_peak_type_kept_when_assay_alone_is_ambiguous(self):
        df = _rows(("ChIP_broad", "ChIP", "H3K4me1", "WT"),
                   ("ChIP_narrow", "ChIP", "H3K4me1", "WT"))
        assert list(plot_levels_labels(df)) == ["WT_ChIP_broad", "WT_ChIP_narrow"]

    def test_mc_collides_on_env_not_assay(self):
        # mC labels carry the context, not the assay, so WGBS and PBAT clash.
        df = _rows(("WGBS", "mC", "WGBS", "Col0"),
                   ("PBAT", "mC", "PBAT", "Col0"),
                   ("WGBS", "mC", "WGBS", "met1"))
        assert list(plot_levels_labels(df)) == ["Col0_WGBS", "Col0_PBAT", "met1"]

    def test_empty(self):
        assert plot_levels_labels(_rows()).empty

    @pytest.mark.parametrize("sheet", [
        "test_samples_CUT.tsv", "test_samples_mC.tsv", "test_samples_colcen.tsv",
        "test_samples_pombe.tsv", "test_samples_hg38_chr21.tsv",
    ])
    def test_real_sheets_have_unique_label_keys(self, sheet):
        # Every (genome, disambiguated levels, mark) a plot label is built from
        # must be unique per analysis sample.
        _, analysis = _analysis(sheet)
        analysis = analysis.copy()
        analysis["levels_label"] = plot_levels_labels(analysis)
        core = analysis.apply(lambda r: r["env"] if r["env"] in ("mC", "sRNA") else r["sample_type"], axis=1)
        keys = list(zip(analysis["ref_genome"], analysis["levels_label"], core))
        assert len(keys) == len(set(keys))

    def test_colcen_pbat_and_wgbs_separated(self):
        _, analysis = _analysis("test_samples_colcen.tsv")
        mc = analysis[(analysis["env"] == "mC") & (analysis["levels_label"] == "Col0_seedling")]
        assert sorted(plot_levels_labels(mc)) == ["Col0_seedling_PBAT", "Col0_seedling_WGBS"]


class TestCheckUniqueLabels:
    def test_unique_passes(self):
        check_unique_labels(["WT_H3K9me2", "WT_H3K9me"], "test")

    def test_duplicate_names_each_once(self):
        with pytest.raises(ValueError, match=r"in the ChIP plots: WT_X\. "):
            check_unique_labels(["WT_X", "WT_Y", "WT_X", "WT_X"], "the ChIP plots")
