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


def _profile(path, groups):
    ticks = ["-2.0Kb", "", "start", "end", "", "2.0Kb"]
    lines = ["bin labels\t\t" + "\t".join(ticks * len(LABELS)),
             "bins\t\t" + "\t".join(str(float(i + 1)) for i in range(NBINS))]
    for label in LABELS:
        for g in groups:
            vals = ["nan" if i == 2 and g.endswith("nostrand") else str(10.0 + i) for i in range(NBINS)]
            lines.append(f"{label}\t{g}\t" + "\t".join(vals))
    path.write_text("\n".join(lines) + "\n")


def _samples(path, labels=LABELS):
    path.write_text("".join(f"{l}\t{l.split('_')[-1]}\t#1f77b4\n" for l in labels))


def _run(tmp_path, groups, labels=LABELS):
    data, samples, out = tmp_path / "p.tab", tmp_path / "s.tsv", tmp_path / "out.pdf"
    _profile(data, groups)
    _samples(samples, labels)
    proc = subprocess.run(["Rscript", str(SCRIPT), str(data), str(samples), str(out), "mean methylation"],
                          capture_output=True, text=True)
    return proc, out


@pytest.mark.parametrize("groups", [["all_genes"], ["all_TEs", "all_TEs_nostrand"]])
def test_draws_a_pdf(tmp_path, groups):
    proc, out = _run(tmp_path, groups)
    assert proc.returncode == 0, proc.stderr
    assert out.read_bytes().startswith(b"%PDF")


def test_column_without_a_context_is_an_error(tmp_path):
    proc, _ = _run(tmp_path, ["all_genes"], labels=LABELS[:-1])
    assert proc.returncode != 0
    assert "ddm1_mCHH" in proc.stderr
