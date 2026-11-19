"""plotHeatmap memory: the lowmem wrapper and the matrix-sized mem_mb request.

matplotlib rasterizes each heatmap of a PDF on a full-figure canvas and the
PDF backend kept a view of every one until the file was written, so
plotHeatmap's memory grew with the square of the column count
(plotheatmap_lowmem.py). The rules run plotHeatmap through the wrapper and
ask for memory sized from the matrix header (heatmap_memory.py).
"""

import gzip
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"
PROFILE = REPO / "profiles" / "default" / "config.yaml"
sys.path.insert(0, str(SCRIPTS))

from heatmap_memory import (  # noqa: E402
    FALLBACK_MB, MIN_MB, heatmap_mem_mb, matrix_dimensions,
)

HEATMAP_RULES = ("plotting_heatmap_on_targetfile", "plotting_sorted_heatmap_on_targetfile")


def header_only(path, regions, columns, bins=180, groups=1):
    """A matrix file with just the header line: all heatmap_mem_mb reads."""
    gb = [round(i * regions / groups) for i in range(groups + 1)]
    header = {"group_labels": [f"g{i}" for i in range(groups)], "group_boundaries": gb,
              "sample_labels": [f"s{i}__mark" for i in range(columns)],
              "sample_boundaries": [i * bins for i in range(columns + 1)]}
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header) + "\n")
    return path


def test_dimensions_from_header(tmp_path):
    path = header_only(tmp_path / "m.gz", 1234, 7, bins=90, groups=2)
    assert matrix_dimensions(path) == (1234, 630, 7)


# Peak RSS measured with deeptools 4.0.0 through the wrapper (180 bins per
# column): the request must cover each one.
MEASURED_GB = [(2000, 60, 0.77), (2000, 120, 1.40), (20000, 30, 3.45), (50000, 60, 16.75)]


@pytest.mark.parametrize("regions,columns,peak_gb", MEASURED_GB)
def test_request_covers_measured_peaks(tmp_path, regions, columns, peak_gb):
    mb = heatmap_mem_mb(header_only(tmp_path / "m.gz", regions, columns))
    assert mb >= peak_gb * 1000 * 1.2
    # and does not ask for far more than the run needs (small runs: the base)
    assert mb <= max(2 * MIN_MB, peak_gb * 1000 * 2.5)


def test_small_matrix_gets_the_floor(tmp_path):
    assert heatmap_mem_mb(header_only(tmp_path / "m.gz", 100, 3)) == MIN_MB


def test_rounded_to_whole_gigabytes(tmp_path):
    assert heatmap_mem_mb(header_only(tmp_path / "m.gz", 30000, 45)) % 1000 == 0


def test_retry_asks_for_more(tmp_path):
    path = header_only(tmp_path / "m.gz", 20000, 60)
    assert heatmap_mem_mb(path, attempt=2) >= 2 * heatmap_mem_mb(path, attempt=1) - 1000


def test_missing_matrix_falls_back(tmp_path):
    """At DAG time the matrix does not exist yet; Snakemake re-evaluates the
    resource before it schedules the job."""
    assert heatmap_mem_mb(tmp_path / "absent.gz") == FALLBACK_MB
    assert heatmap_mem_mb(tmp_path / "absent.gz", attempt=2) == 2 * FALLBACK_MB


def test_unreadable_matrix_falls_back(tmp_path):
    bad = tmp_path / "bad.gz"
    bad.write_text("not gzip")
    assert heatmap_mem_mb(bad) == FALLBACK_MB


def _rule_block(name):
    text = SMK.read_text()
    block = text[text.index(f"rule {name}:"):]
    return block[:block.index("\nrule ", 1)]


@pytest.mark.parametrize("rule", HEATMAP_RULES)
def test_rules_use_wrapper_and_sized_memory(rule):
    block = _rule_block(rule)
    assert "plotHeatmap -m" not in block
    assert re.search(r'python3 "\{params\.plot_heatmap\}" -m \{input\.matrix\}', block)
    assert "plotheatmap_lowmem.py" in block
    assert "mem_mb = lambda wildcards, input, attempt: heatmap_mem_mb(input.matrix, attempt)" in block


@pytest.mark.parametrize("rule", HEATMAP_RULES)
def test_profile_leaves_memory_to_the_rule(rule):
    # set-resources in the workflow profile overrides the rule's resources, so
    # its entry may set runtime but not mem_mb (nor a tier, which carries one).
    entries = [ln.split("#")[0] for ln in PROFILE.read_text().splitlines()
               if ln.strip().startswith(f"{rule}:")]
    resources = [e for e in entries if e.split(":", 1)[1].strip() not in ("1", "")]
    assert len(resources) == 1, entries
    assert "runtime" in resources[0]
    assert "mem_mb" not in resources[0] and "*" not in resources[0]


# --- the wrapper itself, where deeptools is installed --------------------------

def _deeptools_python():
    exe = shutil.which("plotHeatmap")
    if exe is None:
        return None
    py = Path(exe).parent / "python"
    return py if py.exists() else None


DT_PY = _deeptools_python()
needs_deeptools = pytest.mark.skipif(DT_PY is None, reason="deeptools not on PATH")


def _write_matrix(path, regions, columns, bins=36):
    import random
    rnd = random.Random(1)
    header = {"upstream": [0] * columns, "downstream": [0] * columns,
              "body": [bins * 50] * columns, "bin size": [50] * columns,
              "ref point": [None] * columns, "verbose": False, "bin avg type": "mean",
              "missing data as zero": True, "min threshold": None, "max threshold": None,
              "scale": 1, "skip zeros": False, "nan after end": False, "proc number": 1,
              "sort regions": "keep", "sort using": "mean",
              "unscaled 5 prime": [0] * columns, "unscaled 3 prime": [0] * columns,
              "group_labels": ["genes"], "group_boundaries": [0, regions],
              "sample_labels": [f"s{i}__m" for i in range(columns)],
              "sample_boundaries": [i * bins for i in range(columns + 1)]}
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header) + "\n")
        for r in range(regions):
            vals = "\t".join(f"{rnd.random() * 3:.3f}" for _ in range(columns * bins))
            fh.write(f"chr1\t{r * 1000}\t{r * 1000 + 500}\tr{r}\t0\t+\t{vals}\n")


def _images(pdf):
    data = pdf.read_bytes()
    out = []
    for m in re.finditer(rb"/Subtype /Image.*?stream\r?\n", data, re.S):
        end = data.index(b"endstream", m.end())
        out.append(data[m.end():end])
    return out


def _peak_kb(cmd):
    time = shutil.which("time") or "/usr/bin/time"
    if not Path(time).exists():
        subprocess.run(cmd, check=True, capture_output=True)
        return None
    proc = subprocess.run([time, "-v", *cmd], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-2000:]
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", proc.stderr)
    return int(m.group(1)) if m else None


@needs_deeptools
def test_wrapper_output_is_identical_and_smaller(tmp_path):
    columns = 40
    matrix = tmp_path / "m.gz"
    _write_matrix(matrix, 300, columns)
    args = ["-m", str(matrix), "--colorMap", "seismic", "--interpolationMethod", "bilinear",
            "--sortRegions", "descend", "--sortUsing", "mean", "--dpi", "100"]
    plain, low = tmp_path / "plain.pdf", tmp_path / "low.pdf"
    kb_plain = _peak_kb([shutil.which("plotHeatmap"), *args, "-out", str(plain)])
    kb_low = _peak_kb([str(DT_PY), str(SCRIPTS / "plotheatmap_lowmem.py"), *args, "-out", str(low)])
    assert _images(plain) and _images(plain) == _images(low)
    if kb_plain and kb_low:
        assert kb_low < kb_plain / 2, (kb_plain, kb_low)
