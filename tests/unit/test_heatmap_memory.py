"""Heatmap memory: the matrix-sized mem_mb request (heatmap_memory.py).

The heatmaps are drawn in R (R_heatmap.R), one sample column at a time from a
matrix read once; the rule asks for memory sized from the matrix header. The
peaks below were measured on synthetic matrices (180 bins per column) with
/usr/bin/time -v; the R script's own behaviour is tested in test_heatmap.py.
"""

import gzip
import json
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

HEATMAP_RULES = ("plotting_heatmap_on_targetfile",)


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


# Peak RSS of R_heatmap.R (180 bins per column, theme defaults): the request
# must cover each one.
MEASURED_GB = [(2000, 60, 0.92), (2000, 120, 1.42), (20000, 30, 2.56), (50000, 60, 10.52)]


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
def test_rules_draw_in_r_with_sized_memory(rule):
    block = _rule_block(rule)
    assert "plotHeatmap" not in block and "plotheatmap_lowmem.py" not in block
    assert "R_heatmap.R" in block
    assert "mem_mb = lambda wildcards, input, attempt: heatmap_mem_mb(input.matrix, attempt)" in block


def test_wrapper_is_gone():
    assert not (SCRIPTS / "plotheatmap_lowmem.py").exists()
    assert "plotheatmap_lowmem" not in SMK.read_text()


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
