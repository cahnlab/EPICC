"""R_mc_profile_by_context.R: the mC per-group metaplot, faceted on context.

The input mimics `plotProfile --outFileNameData`: two header rows whose
bin-label row repeats once per sample, then one row per sample and region
group. Skipped where Rscript with ggplot2 is not installed.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "workflow" / "scripts" / "R_mc_profile_by_context.R"


def _have_ggplot():
    if not shutil.which("Rscript"):
        return False
    proc = subprocess.run(["Rscript", "-e", "library(ggplot2)"], capture_output=True)
    return proc.returncode == 0


pytestmark = pytest.mark.skipif(not _have_ggplot(), reason="Rscript with ggplot2 not found")

LABELS = ["WT_mCG", "ddm1_mCG", "WT_mCHG", "ddm1_mCHG", "WT_mCHH", "ddm1_mCHH"]
NBINS = 6


def _profile(path, groups, labels=LABELS):
    ticks = ["-2.0Kb", "", "start", "end", "", "2.0Kb"]
    lines = ["bin labels\t\t" + "\t".join(ticks * len(labels)),
             "bins\t\t" + "\t".join(str(float(i + 1)) for i in range(NBINS))]
    for label in labels:
        for g in groups:
            vals = ["nan" if i == 2 and g.endswith("nostrand") else str(10.0 + i) for i in range(NBINS)]
            lines.append(f"{label}\t{g}\t" + "\t".join(vals))
    path.write_text("\n".join(lines) + "\n")


def _samples(path, labels=LABELS, contexts=None, series=None):
    contexts = contexts or [l.split("_")[-1] for l in labels]
    series = series or [l.split("_")[0] for l in labels]
    path.write_text("".join(f"{l}\t{c}\t#1f77b4\t{s}\n" for l, c, s in zip(labels, contexts, series)))


def _run(tmp_path, groups, labels=LABELS, **samples_kw):
    data, samples, out = tmp_path / "p.tab", tmp_path / "s.tsv", tmp_path / "out.pdf"
    _profile(data, groups, labels)
    _samples(samples, labels, **samples_kw)
    proc = subprocess.run(["Rscript", str(SCRIPT), str(data), str(samples), str(out), "mean methylation"],
                          capture_output=True, text=True)
    return proc, out


@pytest.mark.parametrize("groups", [["all_genes"], ["all_TEs", "all_TEs_nostrand"]])
def test_draws_a_pdf(tmp_path, groups):
    proc, out = _run(tmp_path, groups)
    assert proc.returncode == 0, proc.stderr
    assert out.read_bytes().startswith(b"%PDF")


def test_one_context_whose_labels_dropped_it(tmp_path):
    # methylation_contexts: [CG] -- the context is invariant, so
    # build_plot_labels leaves it off the labels; the series is the label.
    proc, out = _run(tmp_path, ["all_genes"], labels=["WT", "ddm1"],
                     contexts=["mCG", "mCG"], series=["WT", "ddm1"])
    assert proc.returncode == 0, proc.stderr
    assert out.read_bytes().startswith(b"%PDF")


def test_missing_series_is_an_error(tmp_path):
    data, samples, out = tmp_path / "p.tab", tmp_path / "s.tsv", tmp_path / "out.pdf"
    _profile(data, ["all_genes"])
    samples.write_text("".join(f"{l}\t{l.split('_')[-1]}\t#1f77b4\n" for l in LABELS))
    proc = subprocess.run(["Rscript", str(SCRIPT), str(data), str(samples), str(out), "y"],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert "series" in proc.stderr


def test_column_without_a_context_is_an_error(tmp_path):
    data, samples, out = tmp_path / "p.tab", tmp_path / "s.tsv", tmp_path / "out.pdf"
    _profile(data, ["all_genes"])
    _samples(samples, LABELS[:-1])
    proc = subprocess.run(["Rscript", str(SCRIPT), str(data), str(samples), str(out), "y"],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert "ddm1_mCHH" in proc.stderr
