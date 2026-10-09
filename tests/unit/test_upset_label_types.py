"""UpSet sample columns must be assigned to their type by lookup, not by regex.

The UpSet scripts colour intersections and the "exclusive mark" violins by type
(histone mark, TSS group, sRNA size class). They used to find a type's columns
with `grep(type, colnames(mat))`, so mark `H3K9me` also claimed every `H3K9me2`
column. `define_samples_for_upset` now emits the label->type map it already
knows while building the labels, and the scripts read their column sets from it.

The function is exec'd out of the .smk and the R lookup is run out of the
scripts, so these tests check the shipped code rather than a copy.
"""

import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"
SCRIPTS = REPO / "workflow" / "scripts"
UPSET_SCRIPTS = ["R_Upset_plot_peaks.R", "R_Upset_plot_TSS.R",
                 "R_Upset_plot_clusters.R"]

# (sample_name, env, sample_type, levels_label, replicates)
ROWS = [
    ("ChIP_broad__WT__H3K9me__Spombe", "ChIP", "H3K9me", "WT", ["rep1", "rep2"]),
    ("ChIP_broad__WT__H3K9me2__Spombe", "ChIP", "H3K9me2", "WT", ["rep1"]),
    ("ChIP_broad__dcr1__H3K9me2__Spombe", "ChIP", "H3K9me2", "dcr1", ["rep1"]),
    ("ATAC__WT__ATAC__Spombe", "ATAC", "ATAC", "WT", ["rep1"]),
    ("RAMPAGE__WT__RAMPAGE__Spombe", "RNA", "RAMPAGE", "WT", ["rep1", "rep2"]),
    ("RAMPAGE__WT_heat__RAMPAGE__Spombe", "RNA", "RAMPAGE", "WT_heat", ["rep1"]),
    ("sRNA__WT__sRNA__Spombe", "sRNA", "sRNA", "WT", ["rep1", "rep2"]),
    ("sRNA__dcr1__sRNA__Spombe", "sRNA", "sRNA", "dcr1", ["rep1"]),
]


def _load_define_samples_for_upset(allreps):
    text = SMK.read_text()
    m = re.search(r"^def define_samples_for_upset\(.*?(?=^\S)", text, re.S | re.M)
    assert m, "define_samples_for_upset not found"

    analysis_samples = pd.DataFrame(
        [{"sample_name": n, "mapped_name": n, "env": e, "sample_type": t,
          "levels_label": lv, "ref_genome": "Spombe"}
         for n, e, t, lv, _ in ROWS])
    reps = {n: [f"{n}_{r}" for r in rr] for n, _, _, _, rr in ROWS}
    ns = {
        "config": {"srna_heatmap_sizes": [], "upset_allreps": allreps,
                   "srna_min_size": 20, "srna_max_size": 24},
        "analysis_samples": analysis_samples,
        "samples": None,
        "RESULTS_DIR": "results",
        "get_replicate_sample_ids": lambda name, df: reps[name],
        "parse_sample_name": lambda sid: {"replicate": sid.rsplit("_", 1)[1]},
        "get_peaktype_for_env": lambda name, env: "broad",
        "get_sample_info_from_name": lambda name, df, field: "PE",
    }
    exec(compile(m.group(0), str(SMK), "exec"), ns)
    return ns["define_samples_for_upset"]


def _call(env, mode, allreps=False):
    fn = _load_define_samples_for_upset(allreps)
    return fn(SimpleNamespace(ref_genome="Spombe", env=env), mode)


def _map(env, allreps=False):
    raw = _call(env, "label_types", allreps)
    return dict(pair.split("=", 1) for pair in raw.split(",")) if raw else {}


def _labels(env, allreps=False):
    return [p.split(":", 1)[0] for p in _call(env, "pairs", allreps)]


class TestLabelTypeMap:
    def test_chip_default(self):
        assert _map("ChIP") == {
            "WT_H3K9me": "H3K9me",
            "WT_H3K9me2": "H3K9me2",
            "dcr1_H3K9me2": "H3K9me2",
        }

    def test_chip_allreps(self):
        assert _map("ChIP", allreps=True) == {
            "WT_H3K9me_rep1": "H3K9me",
            "WT_H3K9me_rep2": "H3K9me",
            "WT_H3K9me2_rep1": "H3K9me2",
            "dcr1_H3K9me2_rep1": "H3K9me2",
        }

    def test_all_chip_spans_chip_and_atac(self):
        m = _map("all_chip")
        assert m["WT_ATAC"] == "ATAC"
        assert m["WT_H3K9me2"] == "H3K9me2"

    @pytest.mark.parametrize("allreps", [False, True])
    def test_rampage_maps_to_levels(self, allreps):
        m = _map("RAMPAGE", allreps)
        expected = ({"WT_rep1": "WT", "WT_rep2": "WT", "WT_heat_rep1": "WT_heat"}
                    if allreps else {"WT": "WT", "WT_heat": "WT_heat"})
        assert m == expected

    def test_srna_labels_carry_size_class(self):
        m = _map("sRNA")
        types = _call("sRNA", "types").split(":")
        assert types == ["20nt", "21nt", "22nt", "23nt", "24nt", "MIRNA", "Others"]
        for prefix in ("WT_rep1", "WT_rep2", "dcr1_rep1"):
            for t in types:
                assert m[f"{prefix}_{t}"] == t
        assert len(m) == 3 * len(types)

    @pytest.mark.parametrize("env,allreps", [
        ("ChIP", False), ("ChIP", True), ("ATAC", False), ("all_chip", True),
        ("RAMPAGE", False), ("RAMPAGE", True)])
    def test_every_label_is_mapped_to_a_listed_type(self, env, allreps):
        m = _map(env, allreps)
        assert set(m) == set(_labels(env, allreps))
        assert set(m.values()) <= set(_call(env, "types", allreps).split(":"))

    def test_rule_passes_map_to_script(self):
        text = SMK.read_text()
        start = text.index("rule plotting_upset_regions:")
        block = text[start:text.index("\nrule ", start + 1)]
        assert 'define_samples_for_upset(wildcards, "label_types")' in block
        assert ('"{params.types}" "{params.label_types}" "{output.plot}"'
                in block)


def _r_type_cols(script, types, label_types, sampleslist):
    """Run the script's own arg parsing and type_cols lookup in base R."""
    lines = (SCRIPTS / script).read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("types<-"))
    head = lines[start:next(i for i, l in enumerate(lines) if l.startswith("output<-"))]
    tc = next(i for i, l in enumerate(lines) if l.startswith("type_cols<-"))
    code = "\n".join([
        "args<-commandArgs(trailingOnly=TRUE)",
        *head,
        "sampleslist<-unlist(strsplit(args[6], ','))",
        lines[tc], lines[tc + 1],
        "for (t in names(type_cols)) cat(t, ':', paste(type_cols[[t]], collapse=','), '\\n', sep='')",
    ])
    out = subprocess.run(
        ["Rscript", "-e", code, "m", "a", "e", ":".join(types),
         ",".join(f"{k}={v}" for k, v in label_types.items()),
         ",".join(sampleslist)],
        check=True, capture_output=True, text=True).stdout
    return {k: (v.split(",") if v else [])
            for k, v in (l.split(":", 1) for l in out.splitlines())}


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not on PATH")
@pytest.mark.parametrize("script", UPSET_SCRIPTS)
class TestRLookup:
    def test_substring_types_do_not_collide(self, script):
        label_types = {"WT_H3K9me": "H3K9me", "WT_H3K9me2": "H3K9me2",
                       "dcr1_H3K9me2": "H3K9me2"}
        got = _r_type_cols(script, ["H3K9me", "H3K9me2"], label_types,
                           ["WT_H3K9me2", "WT_H3K9me", "dcr1_H3K9me2"])
        assert got == {"H3K9me": ["WT_H3K9me"],
                       "H3K9me2": ["WT_H3K9me2", "dcr1_H3K9me2"]}

    def test_only_present_sample_columns_are_used(self, script):
        # A mapped label with no regions is not a column of mat.
        label_types = {"WT_rep1_21nt": "21nt", "WT_rep1_Others": "Others",
                       "WT_rep1_24nt": "24nt"}
        got = _r_type_cols(script, ["21nt", "24nt", "Others"], label_types,
                           ["WT_rep1_Others", "WT_rep1_21nt"])
        assert got == {"21nt": ["WT_rep1_21nt"], "24nt": [],
                       "Others": ["WT_rep1_Others"]}
