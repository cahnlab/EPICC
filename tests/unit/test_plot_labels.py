"""Plot labels in combined analysis must be unique.

Labels are built from Levels plus the mark, not the assay, so two assays
profiling one mark from the same Levels (CUT_RUN + CUT_TAG H3K27me3, WGBS +
PBAT) would share a label. `disambiguate_labels` appends a discriminator to
the colliding ones only. The builders are exec'd out of the .smk and run on
the real test sheets, so these tests check the shipped code.
"""

import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "workflow" / "scripts"))

from sample_sheet import (  # noqa: E402
    add_compat_columns, build_analysis_name, check_unique_labels,
    disambiguate_labels, get_analysis_samples, get_replicate_sample_ids,
    read_sample_sheet,
)

DATA = REPO / "tests" / "integration" / "data"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"


class TestDisambiguateLabels:
    def test_unique_labels_unchanged(self):
        labels = ["WT_H3K9me2", "dcr1_H3K9me2", "WT_H3K4me3"]
        assert disambiguate_labels(labels, ["ChIP_broad", "ChIP_broad", "ChIP_narrow"],
                                   ["a", "b", "c"]) == labels

    def test_assay_without_peak_type(self):
        got = disambiguate_labels(["WT_H3K27me3", "WT_H3K27me3", "mut_H3K27me3"],
                                  ["CUT_RUN_broad", "CUT_TAG_broad", "CUT_RUN_broad"],
                                  ["a", "b", "c"])
        assert got == ["WT_H3K27me3_CUT_RUN", "WT_H3K27me3_CUT_TAG", "mut_H3K27me3"]

    def test_peak_type_when_base_assay_is_ambiguous(self):
        got = disambiguate_labels(["WT_H3K4me1", "WT_H3K4me1"],
                                  ["ChIP_broad", "ChIP_narrow"], ["a", "b"])
        assert got == ["WT_H3K4me1_ChIP_broad", "WT_H3K4me1_ChIP_narrow"]

    def test_sample_id_for_same_assay(self):
        # Two replicates of one group sharing a Replicate_ID.
        got = disambiguate_labels(["WT_rep1_mCG", "WT_rep1_mCG", "WT_rep2_mCG"],
                                  ["WGBS"] * 3, ["WT_a", "WT_b", "WT_c"])
        assert got == ["WT_rep1_mCG_WT_a", "WT_rep1_mCG_WT_b", "WT_rep2_mCG"]

    def test_only_unresolved_entries_fall_back_to_sample_id(self):
        got = disambiguate_labels(["X", "X", "X"], ["WGBS", "WGBS", "PBAT"],
                                  ["s1", "s2", "s3"])
        assert got == ["X_s1", "X_s2", "X_PBAT"]

    def test_stranded_pair_stays_paired(self):
        got = disambiguate_labels(["WT_RNAseq"] * 4, ["RNAseq", "RNAseq", "RAMPAGE", "RAMPAGE"],
                                  ["a", "a", "b", "b"], ["_plus", "_minus", "_plus", "_minus"])
        assert got == ["WT_RNAseq_RNAseq_plus", "WT_RNAseq_RNAseq_minus",
                       "WT_RNAseq_RAMPAGE_plus", "WT_RNAseq_RAMPAGE_minus"]

    def test_suffix_kept_without_collision(self):
        assert disambiguate_labels(["WT", "WT"], ["A", "A"], ["a", "a"],
                                   ["_plus", "_minus"]) == ["WT_plus", "WT_minus"]

    def test_empty(self):
        assert disambiguate_labels([], [], []) == []


class TestCheckUniqueLabels:
    def test_unique_passes(self):
        check_unique_labels(["WT_H3K9me2", "WT_H3K9me"], "test")

    def test_duplicate_names_each_once(self):
        with pytest.raises(ValueError, match=r"in the ChIP plots: WT_X\. "):
            check_unique_labels(["WT_X", "WT_Y", "WT_X", "WT_X"], "the ChIP plots")


# --- The real builders, on the real test sheets ------------------------------

def _sheets(sheet):
    samples = add_compat_columns(read_sample_sheet(str(DATA / sheet)))
    analysis = add_compat_columns(get_analysis_samples(samples))
    # As in the Snakefile: an analysis row is named by its analysis name.
    analysis["sample_name"] = analysis.apply(build_analysis_name, axis=1)
    analysis["mapped_name"] = analysis["sample_name"]
    return samples, analysis


def _builders(sheet, allreps, disambiguate=True):
    samples, analysis = _sheets(sheet)
    config = yaml.safe_load((REPO / "config" / "epicc-options.yaml").read_text())
    config.update(plot_allreps=allreps, upset_allreps=allreps,
                  srna_heatmap_sizes=[21, 24])
    text = SMK.read_text()
    # The plain-Python helpers ahead of the first rule.
    code = text[:text.index("def define_final_stats_output")]

    def replicate(sid):
        return {"replicate": samples.loc[samples["mapped_name"] == sid, "replicate"].iloc[0]}

    ns = {
        "os": os, "pd": pd, "defaultdict": defaultdict,
        "config": config, "samples": samples, "analysis_samples": analysis,
        "RESULTS_DIR": "results", "GENOMES_DIR": "genomes", "REPO_FOLDER": str(REPO),
        "get_replicate_sample_ids": get_replicate_sample_ids,
        "parse_sample_name": replicate,
        "get_peaktype_for_env": lambda name, env: "broad",
        "get_sample_info_from_name": lambda name, df, field: "PE",
        "get_methylation_contexts": lambda: ["CG", "CHG", "CHH"],
        "disambiguate_labels": disambiguate_labels if disambiguate else
            lambda labels, assays, sids, suffixes=None:
                [l + s for l, s in zip(labels, suffixes or [""] * len(labels))],
        "check_unique_labels": check_unique_labels if disambiguate else lambda *a: None,
    }
    exec(compile(code, str(SMK), "exec"), ns)
    return ns, samples["ref_genome"].iloc[0]


SHEETS = ["test_samples_CUT.tsv", "test_samples_mC.tsv", "test_samples_colcen.tsv"]
CASES = (
    [("define_key_for_plots", env, "labels", strand)
     for env in ["all", "most", "ChIP", "ATAC", "RNA", "sRNA", "mC"]
     for strand in ["unstranded", "plus", "nostrand"]]
    + [("define_samples_for_upset", env, "pairs", None)
       for env in ["ChIP", "ATAC", "all_chip", "RAMPAGE", "sRNA"]]
    + [("define_input_for_pca", env, "labels", None)
       for env in ["mCG", "ChIP", "ATAC", "all_chip"]]
)


def _labels(sheet, allreps, fn, env, mode, strand, disambiguate=True):
    ns, genome = _builders(sheet, allreps, disambiguate)
    wc = SimpleNamespace(ref_genome=genome, env=env)
    if strand:
        wc.strand = strand
    out = ns[fn](wc, mode)
    return [p.split(":", 1)[0] for p in out] if mode == "pairs" else out


@pytest.mark.parametrize("allreps", [False, True])
@pytest.mark.parametrize("sheet", SHEETS)
@pytest.mark.parametrize("fn,env,mode,strand", CASES)
def test_builders_make_unique_labels(sheet, allreps, fn, env, mode, strand):
    raw = _labels(sheet, allreps, fn, env, mode, strand, disambiguate=False)
    final = _labels(sheet, allreps, fn, env, mode, strand)
    assert len(final) == len(raw)
    assert len(set(final)) == len(final)
    counts = Counter(raw)
    for r, f in zip(raw, final):
        if counts[r] == 1:
            assert f == r
        else:
            assert f.startswith(r.removesuffix("_plus").removesuffix("_minus"))


def test_cut_run_and_cut_tag_get_their_assay():
    got = _labels("test_samples_CUT.tsv", False, "define_key_for_plots", "ChIP", "labels", "unstranded")
    assert sorted(got) == ["WT_H3K27me3_CUT_RUN", "WT_H3K27me3_CUT_TAG",
                           "WT_TF1_CUT_RUN", "WT_TF1_CUT_TAG"]


def test_colcen_pbat_and_wgbs_separated():
    got = _labels("test_samples_colcen.tsv", False, "define_key_for_plots", "mC", "labels", "unstranded")
    assert "Col0_seedling_mCG_PBAT" in got and "Col0_seedling_mCG_WGBS" in got
    assert "Col0_flower_mCG" in got  # no collision, unchanged


def test_colcen_pca_separates_the_two_inputs_by_sample_id():
    got = _labels("test_samples_colcen.tsv", False, "define_input_for_pca", "ChIP", "labels", None)
    assert "Input_rdr126ddm1_leaf_ChIP_broad_rep1_rdr126ddm1_leaf_CenH3Input_rep1" in got
    assert "Input_rdr126ddm1_leaf_ChIP_broad_rep1_rdr126ddm1_leaf_H3K9me2Input_rep1" in got
    assert "CenH3_Col0_leaf_ChIP_broad_rep1" in got


def test_colliding_samples_keep_one_colour():
    ns, genome = _builders("test_samples_colcen.tsv", False)
    wc = SimpleNamespace(ref_genome=genome, env="mC")
    labels = ns["define_key_for_plots"](wc, "labels")
    colors = dict(zip(labels, ns["define_key_for_plots"](wc, "colors")))
    assert colors["Col0_seedling_mCG_PBAT"] == colors["Col0_seedling_mCG_WGBS"]
    assert colors["Col0_seedling_mCG_PBAT"] != colors["Col0_flower_mCG"]


def test_bigwigs_stay_aligned_with_labels():
    ns, genome = _builders("test_samples_colcen.tsv", True)
    wc = SimpleNamespace(ref_genome=genome, env="mC")
    labels = ns["define_key_for_plots"](wc, "labels")
    bigwigs = ns["define_key_for_plots"](wc, "bigwigs")
    for lab, bw in zip(labels, bigwigs):
        if lab.startswith("Col0_seedling_rep1_mCG_"):
            assay = lab.rsplit("_", 1)[1]
            assert bw.endswith(f"Col0_seedling_{assay}_rep1__ColCEN__CG.bw")


def test_srna_label_types_follow_final_labels():
    ns, genome = _builders("test_samples_colcen.tsv", False)
    wc = SimpleNamespace(ref_genome=genome, env="sRNA")
    prefixes = [p.split(":", 1)[0] for p in ns["define_samples_for_upset"](wc, "pairs")]
    pairs = ns["define_samples_for_upset"](wc, "label_types").split(",")
    assert {p.rsplit("=", 1)[0].rsplit("_", 1)[0] for p in pairs} == set(prefixes)
