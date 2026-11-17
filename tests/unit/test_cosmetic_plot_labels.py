"""Plot labels are cosmetic: changing them reruns plotting jobs only.

Snakemake reruns a job whenever one of its params changes. Labels used to be
params of the steps that compute data -- computeMatrix --samplesLabel,
multiBigwigSummary/multiBamSummary -l, the UpSet peak merge -- so any change
to how labels are spelled rebuilt every matrix and summary. Data steps now
name samples by stable keys built from sample ids, and labels and colours are
joined on by the plotting rules alone.

The builders are exec'd out of the .smk and the shell snippets are run out of
the rule bodies, so these tests check the shipped code.
"""

import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "unit"))
sys.path.insert(0, str(REPO / "workflow" / "scripts"))
sys.path.insert(0, str(REPO / "workflow"))  # the .smk imports scripts.*

import test_plot_labels as tpl  # noqa: E402
from matrix_sample_labels import display_rows  # noqa: E402

SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"
LABELLER = REPO / "workflow" / "scripts" / "matrix_sample_labels.py"
TEXT = SMK.read_text()

# Rules that compute data: none of them may take a label or a colour.
DATA_RULES = [
    "making_stranded_matrix_on_targetfile", "merging_matrix",
    "computing_matrix_scales", "sort_heatmap", "combine_clusterfiles",
    "combine_peakfiles", "combine_TSS", "prep_browser_on_region",
    "summarize_tracks_pca",
]
PLOT_RULES = [
    "plotting_heatmap_on_targetfile", "plotting_sorted_heatmap_on_targetfile",
    "plotting_profile_on_targetfile", "make_single_loci_browser_plot",
    "plot_PCA_correlation", "plotting_upset_regions",
]
LABEL_MODES = ['"labels"', '"colors"', '"display"', '"browser_display"',
               '"relabel"', '"label_types"']


def _rule(name):
    start = TEXT.index(f"rule {name}:")
    nxt = re.search(r"^(?:rule|checkpoint) ", TEXT[start + 1:], re.M)
    return TEXT[start:start + 1 + nxt.start()] if nxt else TEXT[start:]


def _shell(name):
    return re.search(r'shell:\s*"""(.*?)"""', _rule(name), re.S).group(1)


@pytest.mark.parametrize("rule", DATA_RULES)
def test_data_rules_take_no_labels(rule):
    block = _rule(rule)
    for mode in LABEL_MODES:
        assert mode not in block, f"{rule} depends on {mode}"


@pytest.mark.parametrize("rule", PLOT_RULES)
def test_plot_rules_apply_labels(rule):
    block = _rule(rule)
    assert any(mode in block for mode in LABEL_MODES)


# --- Keys do not move when labels do ----------------------------------------

def _relabelled_builders(sheet, allreps):
    ns, genome = tpl._builders(sheet, allreps)
    real = ns["disambiguate_labels"]
    ns["disambiguate_labels"] = lambda *a, **k: [x + "_RELABELLED" for x in real(*a, **k)]
    return ns, genome


KEY_CASES = (
    [("define_key_for_plots", env, mode, strand)
     for env in ["all", "most", "mC", "RNA", "sRNA"]
     for strand in ["unstranded", "plus", "nostrand"]
     for mode in ["keys", "bigwigs", "tracks"]]
    + [("define_input_for_pca", env, mode, None)
       for env in ["mCG", "ChIP", "all_chip"] for mode in ["keys", "tracks"]]
    + [("define_samples_for_upset", env, mode, None)
       for env in ["ChIP", "all_chip", "sRNA"] for mode in ["pairs", "files"]]
)


@pytest.mark.parametrize("allreps", [False, True])
@pytest.mark.parametrize("sheet", ["test_samples_colcen.tsv", "test_samples_CUT.tsv"])
@pytest.mark.parametrize("fn,env,mode,strand", KEY_CASES)
def test_data_inputs_ignore_labels(sheet, allreps, fn, env, mode, strand):
    wc = SimpleNamespace(env=env)
    if strand:
        wc.strand = strand
    out = []
    for build in (tpl._builders, _relabelled_builders):
        ns, genome = build(sheet, allreps)
        wc.ref_genome = genome
        out.append(ns[fn](wc, mode))
        if mode in ("keys", "pairs"):
            assert not any("RELABELLED" in x for x in out[-1])
    assert out[0] == out[1]


@pytest.mark.parametrize("env", ["all", "most", "mC", "RNA", "sRNA"])
def test_display_covers_both_strand_layouts(env):
    ns, genome = tpl._builders("test_samples_colcen.tsv", False)
    rows = {}
    for line in ns["define_key_for_plots"](SimpleNamespace(ref_genome=genome, env=env), "display").split("\n"):
        key, label, color, track = line.split("\t")
        rows[key] = (label, color, track)
    for strand in ("unstranded", "plus"):
        wc = SimpleNamespace(ref_genome=genome, env=env, strand=strand)
        fn = ns["define_key_for_plots"]
        for key, label, color, track in zip(fn(wc, "keys"), fn(wc, "labels"), fn(wc, "colors"), fn(wc, "tracks")):
            assert rows[key] == (label, color, track)


def test_keys_end_in_their_track():
    # computing_matrix_scales reads the mark back from the last '__' field, and
    # prep_browser_on_region negates keys ending in '_minus'.
    ns, genome = tpl._builders("test_samples_colcen.tsv", True)
    wc = SimpleNamespace(ref_genome=genome, env="all", strand="unstranded")
    keys, tracks = ns["define_key_for_plots"](wc, "keys"), ns["define_key_for_plots"](wc, "tracks")
    assert any(t.endswith("_minus") for t in tracks)
    for key, track in zip(keys, tracks):
        assert key.rsplit("__", 1)[1] == track


# --- matrix_sample_labels.py ------------------------------------------------

def _write_matrix(path, sample_keys, nregions=4, values=None):
    nbins = 3
    n = len(sample_keys)
    # nbins = 3: 200 bp upstream + 100 bp downstream in 100 bp bins. The
    # per-sample fields carry one entry per sample.
    header = {
        "upstream": [200] * n, "downstream": [100] * n, "body": [0] * n,
        "bin size": [100] * n, "ref point": ["TSS"] * n, "verbose": False,
        "bin avg type": "mean", "missing data as zero": False,
        "min threshold": None, "max threshold": None, "scale": 1,
        "skip zeros": False, "nan after end": False, "proc number": 1,
        "sort regions": "keep", "sort using": "mean",
        "unscaled 5 prime": [0] * n, "unscaled 3 prime": [0] * n,
        "group_labels": ["genes"], "group_boundaries": [0, nregions],
        "sample_labels": list(sample_keys),
        "sample_boundaries": [i * nbins for i in range(n + 1)],
    }
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header) + "\n")
        for r in range(nregions):
            row = ["chr1", str(r * 1000), str(r * 1000 + 500), f"g{r}", "0", "+"]
            for s in range(len(sample_keys)):
                v = values[s] if values else s + 1
                row += [f"{v + 0.1 * (r + b) if v else 0:.2f}" for b in range(nbins)]
            fh.write("\t".join(row) + "\n")


def test_display_rows_follow_matrix_order(tmp_path):
    matrix, display = tmp_path / "m.gz", tmp_path / "d.tsv"
    _write_matrix(matrix, ["s2__mCG", "s1__mCG", "s1__mCHG"])
    display.write_text("s1__mCG\tWT_mCG\t#111111\tmCG\n"
                       "s1__mCHG\tWT_mCHG\t#111111\tmCHG\n"
                       "s2__mCG\tddm1_mCG\t#222222\tmCG\n"
                       "s3__mCG\tunused\t#333333\tmCG\n")
    rows = display_rows(matrix, display)
    assert [r.split("\t")[1] for r in rows] == ["ddm1_mCG", "WT_mCG", "WT_mCHG"]


def test_unknown_column_is_an_error(tmp_path):
    matrix, display = tmp_path / "m.gz", tmp_path / "d.tsv"
    _write_matrix(matrix, ["s1__mCG", "old__mCG"])
    display.write_text("s1__mCG\tWT\t#111111\tmCG\n")
    proc = subprocess.run([sys.executable, str(LABELLER), str(matrix), str(display)],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert "old__mCG" in proc.stderr


# --- Shell snippets out of the plotting rules --------------------------------

def _bash(script, tmp_path):
    tmpdir = tmp_path / "tmpdir"
    tmpdir.mkdir(exist_ok=True)
    return subprocess.run(["bash", "-c", "set -euo pipefail\n" + script], capture_output=True, text=True,
                          env={"PATH": os.environ["PATH"], "TMPDIR": str(tmpdir)})


def _fill(body, values):
    for placeholder, value in values.items():
        body = body.replace(placeholder, value)
    return body.replace("{{", "{").replace("}}", "}")


def test_browser_plot_joins_labels_onto_extracted_tracks(tmp_path):
    body = _shell("make_single_loci_browser_plot")
    snippet = body[body.index("printf '%s\\n' \"{params.display}\""):body.index("Rscript")]
    filenames = tmp_path / "filenames.txt"
    filenames.write_text("Name\tPath\tYmin\tYmax\n"
                         "a__RNAseq_plus\t/t/a__RNAseq_plus.bw\t0\t5\n"
                         "a__RNAseq_minus\t/t/a__RNAseq_minus.bw\t-5\t0\n")
    display = ("a__RNAseq_minus\tWT_RNAseq_minus\t#back\t#track\t#plus\t#minus\n"
               "a__RNAseq_plus\tWT_RNAseq_plus\t#back\t#track\t#plus\t#minus")
    proc = _bash(_fill(snippet, {"{params.display}": display, "{input.filenames}": str(filenames)})
                 + 'cat "$TMPDIR/filenames.tsv"', tmp_path)
    assert proc.returncode == 0, proc.stderr
    lines = [l.split("\t") for l in proc.stdout.splitlines()]
    assert lines[0] == ["Name", "Path", "Backcolor", "Trackcolor", "Fillcolorplus",
                        "Fillcolorminus", "Ymin", "Ymax"]
    assert lines[1] == ["WT_RNAseq_plus", "/t/a__RNAseq_plus.bw", "#back", "#track",
                        "#plus", "#minus", "0", "5"]
    assert lines[2][0] == "WT_RNAseq_minus" and lines[2][6:] == ["-5", "0"]

    filenames.write_text("Name\tPath\tYmin\tYmax\nstale__H3K9me2\t/t/x.bw\t0\t1\n")
    proc = _bash(_fill(snippet, {"{params.display}": display, "{input.filenames}": str(filenames)}), tmp_path)
    assert proc.returncode != 0 and "stale__H3K9me2" in proc.stderr


def test_upset_plot_relabels_merged_samples(tmp_path):
    body = _shell("plotting_upset_regions")
    snippet = body[:body.index("Rscript")]
    merged = tmp_path / "merged.bed"
    merged.write_text("Chr\tStart\tStop\tPeakID\tSamples\n"
                      "chr1\t1\t10\tp1\tA__G,B__G\n"
                      "chr1\t20\t30\tp2\tB__G\n")
    relabel = "A__G\tWT_H3K9me\nB__G\tWT_H3K9me2"
    proc = _bash(_fill(snippet, {"{params.relabel}": relabel, "{input.mergedfile}": str(merged)})
                 + 'cat "$TMPDIR/merged_labelled.bed"', tmp_path)
    assert proc.returncode == 0, proc.stderr
    rows = [l.split("\t") for l in proc.stdout.splitlines()]
    assert rows[0][4] == "Samples"
    assert rows[1][4] == "WT_H3K9me,WT_H3K9me2" and rows[2][4] == "WT_H3K9me2"


@pytest.mark.requires_deeptools
@pytest.mark.skipif(shutil.which("computeMatrixOperations") is None, reason="deeptools not found")
def _scales(tmp_path, mode, keys, values):
    body = _shell("computing_matrix_scales")
    matrix = tmp_path / "final.gz"
    _write_matrix(matrix, keys, values=values)
    out = {k: tmp_path / f"{k}.txt" for k in
           ("params_heatmap", "params_profile", "params_regions", "temp_values",
            "temp_profile_values")}
    script = _fill(body, {
        **{f"{{output.{k}}}": str(v) for k, v in out.items()},
        "{output.temp_profile}": str(tmp_path / "p.pdf"),
        "{input.matrix}": str(matrix), "{log}": str(tmp_path / "log"),
        "{params.labeller}": str(REPO / "workflow" / "scripts" / "matrix_region_labels.py"),
        "{params.scales}": mode, "{params.profile}": "mean",
        "{params.matrix}": "tss", "{params.env}": "ChIP", "{params.target_name}": "genes",
        "{params.ref_genome}": "G", "{params.cg_scale}": "1", "{params.chg_scale}": "1",
        "{params.chh_scale}": "1",
    })
    proc = _bash(script, tmp_path)
    assert proc.returncode == 0, proc.stderr
    heat = out["params_heatmap"].read_text().split()
    return heat[heat.index("--zMin") + 1:heat.index("--zMax")], heat[heat.index("--zMax") + 1:]


@pytest.mark.requires_deeptools
@pytest.mark.skipif(shutil.which("computeMatrixOperations") is None, reason="deeptools not found")
def test_type_scales_group_columns_by_exact_mark(tmp_path):
    # H3K9me must not pick up H3K9me2's columns, as a substring grep would.
    zmin, zmax = _scales(tmp_path, "type", ["a__H3K9me", "b__H3K9me", "a__H3K9me2"], [1, 2, 50])
    assert len(zmin) == len(zmax) == 3
    # Columns 1-2 (H3K9me) share a range that excludes column 3 (H3K9me2, ~50).
    assert zmax[0] == zmax[1] and float(zmax[0]) < 10 < float(zmin[2])


@pytest.mark.requires_deeptools
@pytest.mark.skipif(shutil.which("computeMatrixOperations") is None, reason="deeptools not found")
def test_sample_scales_one_range_per_column(tmp_path):
    # An all-zero mCHG column falls back to the fixed mCHG scale, found by the
    # key's track rather than by the text of a label.
    zmin, zmax = _scales(tmp_path, "sample", ["a__mCG", "a__mCHG", "b__mCG"], [5, 0, 7])
    assert len(zmin) == len(zmax) == 3
    assert zmin[1] == "0" and zmax[1] == "1"
    assert float(zmax[0]) < float(zmax[2])
