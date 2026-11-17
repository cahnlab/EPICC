"""Plot labels: show a component only where it varies, and stay unique.

A label is built from components -- each Levels factor value, the mark
(IP_target / mC context / RNA type), the assay, and the replicate where
replicates are plotted individually. `build_plot_labels` shows a component on
every label if it differs between any of the samples plotted together, and
drops it if it is the same for all of them. Labels still identical after that
get the full assay token, then the Sample_ID. The builders are exec'd out of
the .smk and run on the real test sheets, so these tests check the shipped
code.
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
sys.path.insert(0, str(REPO / "workflow"))  # the .smk imports scripts.*

from sample_sheet import (  # noqa: E402
    add_compat_columns, build_analysis_name, build_plot_labels,
    check_unique_labels, get_analysis_samples, get_replicate_sample_ids,
    plot_label_parts, read_sample_sheet,
)

DATA = REPO / "tests" / "integration" / "data"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"


def E(levels, mark="", assay="ChIP_broad", replicate="", sid=None, suffix=""):
    return dict(levels=levels, mark=mark, assay=assay, replicate=replicate,
                sid=sid or f"{levels}_{mark}_{assay}_{replicate}", suffix=suffix)


class TestBuildPlotLabels:
    def test_varying_components_shown_invariant_dropped(self):
        got = build_plot_labels([
            E("genotype:WT,tissue:leaf", "H3K9me2"),
            E("genotype:ddm1,tissue:leaf", "H3K9me2"),
            E("genotype:WT,tissue:leaf", "H3K27me3"),
        ])
        # tissue and assay are the same everywhere; genotype and mark are not.
        assert got == ["WT_H3K9me2", "ddm1_H3K9me2", "WT_H3K27me3"]

    def test_varying_factor_is_shown_on_every_label(self):
        got = build_plot_labels([
            E("genotype:Col0,tissue:flower", "mCG", "EMseq"),
            E("genotype:Col0,tissue:seedling", "mCG", "WGBS"),
            E("genotype:Col0,tissue:seedling", "mCG", "PBAT"),
        ])
        # Old rule: only the colliding pair got an assay. Now assay varies, so
        # every label carries it; genotype and mark don't, so none does.
        assert got == ["flower_EMseq", "seedling_WGBS", "seedling_PBAT"]

    def test_single_sample_keeps_its_levels(self):
        assert build_plot_labels([E("genotype:WT,tissue:root", "H3K9me2")]) == ["WT_root"]

    def test_all_components_invariant_keeps_levels(self):
        got = build_plot_labels([E("genotype:WT", "H3K9me2", suffix="_plus"),
                                 E("genotype:WT", "H3K9me2", suffix="_minus")])
        assert got == ["WT_plus", "WT_minus"]

    def test_levels_dropped_when_only_the_mark_varies(self):
        got = build_plot_labels([E("genotype:WT", "H3K9me2"), E("genotype:WT", "H3K27me3")])
        assert got == ["H3K9me2", "H3K27me3"]

    def test_assay_without_peak_type(self):
        got = build_plot_labels([
            E("genotype:WT", "H3K27me3", "CUT_RUN_broad"),
            E("genotype:WT", "H3K27me3", "CUT_TAG_broad"),
            E("genotype:mut", "H3K27me3", "CUT_RUN_broad"),
        ])
        assert got == ["WT_CUT_RUN", "WT_CUT_TAG", "mut_CUT_RUN"]

    def test_peak_type_when_only_it_varies(self):
        got = build_plot_labels([E("genotype:WT", "H3K4me1", "ChIP_broad"),
                                 E("genotype:WT", "H3K4me1", "ChIP_narrow")])
        assert got == ["ChIP_broad", "ChIP_narrow"]

    def test_peak_type_breaks_a_tie_the_base_assay_left(self):
        got = build_plot_labels([
            E("genotype:WT", "H3K4me1", "ChIP_broad"),
            E("genotype:WT", "H3K4me1", "ChIP_narrow"),
            E("genotype:WT", "H3K4me1", "CUT_TAG_broad"),
        ])
        assert got == ["ChIP_broad", "ChIP_narrow", "CUT_TAG"]

    def test_assay_not_repeated_when_the_mark_names_it(self):
        got = build_plot_labels([
            E("genotype:WT", "H3K9me2", "ChIP_broad"),
            E("genotype:WT", "ATAC", "ATAC"),
            E("genotype:WT", "RNAseq", "RNAseq"),
            E("genotype:WT", "RAMPAGE", "RAMPAGE"),
            E("genotype:WT", "sRNA_24nt", "sRNA"),
        ])
        assert got == ["H3K9me2_ChIP", "ATAC", "RNAseq", "RAMPAGE", "sRNA_24nt"]

    def test_replicate_shown_when_it_varies(self):
        got = build_plot_labels([E("genotype:WT", "mCG", "WGBS", "rep1"),
                                 E("genotype:WT", "mCG", "WGBS", "rep2"),
                                 E("genotype:ddm1", "mCG", "WGBS", "rep1")])
        assert got == ["WT_rep1", "WT_rep2", "ddm1_rep1"]

    def test_replicate_dropped_when_every_sample_has_the_same(self):
        got = build_plot_labels([E("genotype:WT", "mCG", "WGBS", "rep1"),
                                 E("genotype:ddm1", "mCG", "WGBS", "rep1")])
        assert got == ["WT", "ddm1"]

    def test_sample_id_last_resort(self):
        # Two Inputs from one Levels serving different marks.
        got = build_plot_labels([
            E("genotype:WT", "Input", replicate="rep1", sid="WT_ChIP_Input_rep1"),
            E("genotype:WT", "Input", replicate="rep1", sid="WT_ChIP_Input_rep1b"),
            E("genotype:WT", "H3K9me2", replicate="rep1"),
        ])
        assert got == ["Input_ChIP_broad_WT_ChIP_Input_rep1",
                       "Input_ChIP_broad_WT_ChIP_Input_rep1b", "H3K9me2"]

    def test_tiebreak_on_an_all_invariant_label_keeps_levels(self):
        got = build_plot_labels([E("genotype:WT", "Input", sid="a"),
                                 E("genotype:WT", "Input", sid="b")])
        assert got == ["WT_ChIP_broad_a", "WT_ChIP_broad_b"]

    def test_stranded_pair_stays_paired(self):
        got = build_plot_labels([
            E("genotype:WT", "RNAseq", "RNAseq", sid="a", suffix="_plus"),
            E("genotype:WT", "RNAseq", "RNAseq", sid="a", suffix="_minus"),
            E("genotype:dcr1", "RNAseq", "RNAseq", sid="b", suffix="_plus"),
            E("genotype:dcr1", "RNAseq", "RNAseq", sid="b", suffix="_minus"),
        ])
        assert got == ["WT_plus", "WT_minus", "dcr1_plus", "dcr1_minus"]

    def test_omit_mark_names_a_line_across_context_panels(self):
        entries = [E("genotype:WT", "mCG", "WGBS"), E("genotype:ddm1", "mCG", "WGBS"),
                   E("genotype:WT", "mCHG", "WGBS"), E("genotype:ddm1", "mCHG", "WGBS")]
        assert build_plot_labels(entries) == ["WT_mCG", "ddm1_mCG", "WT_mCHG", "ddm1_mCHG"]
        assert build_plot_labels(entries, omit=("mark",)) == ["WT", "ddm1", "WT", "ddm1"]

    def test_omit_never_leaves_an_empty_name(self):
        entries = [E("genotype:WT", "mCG", "WGBS"), E("genotype:WT", "mCHG", "WGBS")]
        assert build_plot_labels(entries) == ["mCG", "mCHG"]
        assert build_plot_labels(entries, omit=("mark",)) == ["WT", "WT"]

    def test_parts_are_structured(self):
        parts = plot_label_parts([E("genotype:WT", "mCG", "WGBS", "rep1"),
                                  E("genotype:ddm1", "mCG", "PBAT", "rep2")])
        assert parts[0] == [("level", "WT"), ("assay", "WGBS"), ("replicate", "rep1")]

    def test_empty(self):
        assert build_plot_labels([]) == []


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


def _builders(sheet, allreps, contexts=("CG", "CHG", "CHH")):
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
        "get_methylation_contexts": lambda: list(contexts),
        "build_plot_labels": build_plot_labels,
        "check_unique_labels": check_unique_labels,
    }
    exec(compile(code, str(SMK), "exec"), ns)
    return ns, samples["ref_genome"].iloc[0]


SHEETS = ["test_samples_CUT.tsv", "test_samples_mC.tsv", "test_samples_colcen.tsv"]
CASES = (
    [("define_key_for_plots", env, "labels", strand)
     for env in ["all", "most", "ChIP", "ATAC", "RNA", "sRNA", "mC"]
     for strand in ["unstranded", "plus", "nostrand"]]
    + [("define_samples_for_upset", env, "labels", None)
       for env in ["ChIP", "ATAC", "all_chip", "RAMPAGE", "sRNA"]]
    + [("define_input_for_pca", env, "labels", None)
       for env in ["mCG", "ChIP", "ATAC", "all_chip"]]
)


def _labels(sheet, allreps, fn, env, mode, strand, contexts=("CG", "CHG", "CHH")):
    ns, genome = _builders(sheet, allreps, contexts)
    wc = SimpleNamespace(ref_genome=genome, env=env)
    if strand:
        wc.strand = strand
    return ns[fn](wc, mode)


@pytest.mark.parametrize("allreps", [False, True])
@pytest.mark.parametrize("sheet", SHEETS)
@pytest.mark.parametrize("fn,env,mode,strand", CASES)
def test_builders_make_unique_nonempty_labels(sheet, allreps, fn, env, mode, strand):
    labels = _labels(sheet, allreps, fn, env, mode, strand)
    assert len(set(labels)) == len(labels)
    assert all(lab and not lab.startswith("_") and "__" not in lab for lab in labels)


def test_cut_run_and_cut_tag_get_their_assay():
    got = _labels("test_samples_CUT.tsv", False, "define_key_for_plots", "ChIP", "labels", "unstranded")
    # Every sample is WT, so the Levels are dropped; mark and assay both vary.
    assert sorted(got) == ["H3K27me3_CUT_RUN", "H3K27me3_CUT_TAG", "TF1_CUT_RUN", "TF1_CUT_TAG"]


def test_colcen_mc_labels_show_every_varying_component():
    got = _labels("test_samples_colcen.tsv", False, "define_key_for_plots", "mC", "labels", "unstranded")
    # Assays vary across the mC samples, so every label carries one -- not just
    # the PBAT/WGBS pair that would otherwise collide.
    assert "Col0_seedling_mCG_PBAT" in got and "Col0_seedling_mCG_WGBS" in got
    assert "Col0_flower_mCG_EMseq" in got


def test_colcen_rna_drops_invariant_tissue_and_mark():
    got = _labels("test_samples_colcen.tsv", False, "define_key_for_plots", "RNA", "labels", "plus")
    assert got == ["Col0", "suvh13"]


def test_colcen_cg_only_drops_the_context():
    got = _labels("test_samples_colcen.tsv", False, "define_key_for_plots", "mC", "labels",
                  "unstranded", contexts=("CG",))
    assert "Col0_flower_EMseq" in got and not any("mCG" in lab for lab in got)


def test_colcen_pca_inputs_are_distinct_without_a_suffix():
    # The CenH3 and H3K9me2 Inputs come from different genotypes (rdr126 vs
    # rdr126ddm1, per GSE132005), so their labels differ on Levels alone.
    got = _labels("test_samples_colcen.tsv", False, "define_input_for_pca", "ChIP", "labels", None)
    assert "rdr126_leaf_Input_ChIP_rep1" in got
    assert "rdr126ddm1_leaf_Input_ChIP_rep1" in got
    assert "Col0_leaf_CenH3_ChIP_rep1" in got


def test_colours_follow_the_full_levels_label():
    # Colours are keyed on levels_label, not the shortened label: samples of one
    # Levels share a colour whatever their labels show.
    ns, genome = _builders("test_samples_colcen.tsv", False)
    wc = SimpleNamespace(ref_genome=genome, env="mC")
    labels = ns["define_key_for_plots"](wc, "labels")
    colors = dict(zip(labels, ns["define_key_for_plots"](wc, "colors")))
    assert colors["Col0_seedling_mCG_PBAT"] == colors["Col0_seedling_mCG_WGBS"]
    assert colors["Col0_seedling_mCG_PBAT"] != colors["Col0_flower_mCG_EMseq"]


def test_colour_is_keyed_on_levels_label_not_the_label():
    # RNA labels drop the invariant tissue ('Col0'), but the colour still comes
    # from the full levels_label ('Col0_seedling').
    ns, genome = _builders("test_samples_colcen.tsv", False)
    key = ns["define_key_for_plots"]
    wc = SimpleNamespace(ref_genome=genome, env="RNA", strand="plus")
    colors = dict(zip(key(wc, "labels"), key(wc, "colors")))
    palette = ns["assign_colors"](["Col0_seedling", "suvh13_seedling"])
    assert colors == {"Col0": palette["Col0_seedling"], "suvh13": palette["suvh13_seedling"]}


def test_bigwigs_stay_aligned_with_labels():
    ns, genome = _builders("test_samples_colcen.tsv", True)
    wc = SimpleNamespace(ref_genome=genome, env="mC")
    labels = ns["define_key_for_plots"](wc, "labels")
    bigwigs = ns["define_key_for_plots"](wc, "bigwigs")
    seen = 0
    for lab, bw in zip(labels, bigwigs):
        if lab.startswith("Col0_seedling_mCG_") and lab.endswith("_rep1"):
            assay = lab.split("_")[3]
            assert bw.endswith(f"Col0_seedling_{assay}_rep1__ColCEN__CG.bw")
            seen += 1
    assert seen == 2  # PBAT and WGBS


def test_srna_label_types_follow_final_labels():
    ns, genome = _builders("test_samples_colcen.tsv", False)
    wc = SimpleNamespace(ref_genome=genome, env="sRNA")
    prefixes = ns["define_samples_for_upset"](wc, "labels")
    pairs = ns["define_samples_for_upset"](wc, "label_types").split(",")
    assert {p.rsplit("=", 1)[0].rsplit("_", 1)[0] for p in pairs} == set(prefixes)


# --- The context of each mC matrix column -----------------------------------

def test_tracks_name_each_columns_context():
    # R_mc_profile_by_context.R facets the mC per-group profile on these.
    ns, genome = _builders("test_samples_colcen.tsv", True)
    wc = SimpleNamespace(ref_genome=genome, env="mC")
    key = ns["define_key_for_plots"]
    labels, tracks, bigwigs = key(wc, "labels"), key(wc, "tracks"), key(wc, "bigwigs")
    assert len(tracks) == len(labels) == len(bigwigs)
    assert set(tracks) == {"mCG", "mCHG", "mCHH"}
    for lab, ctx, bw in zip(labels, tracks, bigwigs):
        assert f"_{ctx}" in lab and bw.endswith(f"__{ctx[1:]}.bw")


def test_series_name_each_sample_across_context_panels():
    ns, genome = _builders("test_samples_colcen.tsv", False)
    wc = SimpleNamespace(ref_genome=genome, env="mC")
    key = ns["define_key_for_plots"]
    labels, series, tracks = key(wc, "labels"), key(wc, "series"), key(wc, "tracks")
    for lab, ser, ctx in zip(labels, series, tracks):
        assert lab.replace(f"_{ctx}", "") == ser
    # One line per sample: each series appears once per context.
    assert set(Counter(series).values()) == {3}


def test_series_with_one_context_is_the_label():
    ns, genome = _builders("test_samples_colcen.tsv", True, contexts=("CG",))
    wc = SimpleNamespace(ref_genome=genome, env="mC")
    key = ns["define_key_for_plots"]
    assert key(wc, "series") == key(wc, "labels")
    assert not any("mCG" in lab for lab in key(wc, "labels"))
