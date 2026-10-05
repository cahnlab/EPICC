"""Tests for workflow/scripts/merge_cx_reports.sh.

Replicate reports that line up row for row are summed side by side; reports
that don't are merged by sort + bedtools merge. Both must give the same
counts.
"""

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "scripts" / "merge_cx_reports.sh"

if not (shutil.which("bedtools") and shutil.which("pigz")):
    pytest.skip("bedtools/pigz not found", allow_module_level=True)

# Genome order, not sorted order: chr2 before chr10 before chr1.
SITES = [
    ("chr2", 3, "+", "CG", "CGA"),
    ("chr2", 4, "-", "CG", "CGT"),
    ("chr2", 40, "+", "CHH", "CTA"),
    ("chr10", 7, "-", "CHG", "CAG"),
    ("chr1", 5, "+", "CHH", "CCT"),
    ("chr1", 9, "+", "CG", "CGG"),
]


def _report(path, counts, sites=SITES):
    with gzip.open(path, "wt") as fh:
        for (chrom, pos, strand, ctx, tri), (m, u) in zip(sites, counts):
            fh.write(f"{chrom}\t{pos}\t{strand}\t{m}\t{u}\t{ctx}\t{tri}\n")
    return str(path)


def _merge(tmp_path, reports):
    out = tmp_path / "merged.CX_report.txt.gz"
    result = subprocess.run(
        ["bash", str(SCRIPT), "2", str(out), *reports],
        capture_output=True, text=True, env={**os.environ, "TMPDIR": str(tmp_path)},
    )
    assert result.returncode == 0, result.stderr
    with gzip.open(out, "rt") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return rows, result.stdout


def _by_site(rows):
    return {(r[0], int(r[1])): r for r in rows}


def test_aligned_replicates_are_summed_in_order(tmp_path):
    a = _report(tmp_path / "a.gz", [(1, 0), (0, 0), (2, 3), (0, 4), (5, 5), (1, 1)])
    b = _report(tmp_path / "b.gz", [(2, 1), (0, 0), (1, 1), (3, 0), (0, 2), (0, 0)])
    rows, stdout = _merge(tmp_path, [a, b])
    assert "merging by sort" not in stdout
    assert [r[:2] for r in rows] == [[c, str(p)] for c, p, *_ in SITES]
    assert rows[0] == ["chr2", "3", "+", "3", "1", "CG", "CGA"]
    assert rows[2] == ["chr2", "40", "+", "3", "4", "CHH", "CTA"]


def test_chromosome_order_differing_between_replicates(tmp_path):
    # Bismark 3.x writes chromosomes in a different order on every run.
    counts = [(1, 0), (0, 0), (2, 3), (0, 4), (5, 5), (1, 1)]
    a = _report(tmp_path / "a.gz", counts)
    order = [4, 5, 3, 0, 1, 2]  # chr1, chr10, chr2
    b = _report(tmp_path / "b.gz", [counts[i] for i in order],
                [SITES[i] for i in order])
    rows, stdout = _merge(tmp_path, [a, b])
    assert "merging by sort" not in stdout
    assert _by_site(rows) == {
        (c, p): [c, str(p), s, str(2 * m), str(2 * u), ctx, tri]
        for (c, p, s, ctx, tri), (m, u) in zip(SITES, counts)
    }


def test_aligned_matches_sort_merge(tmp_path):
    counts = [[(1, 0), (0, 0), (2, 3), (0, 4), (5, 5), (1, 1)],
              [(2, 1), (0, 0), (1, 1), (3, 0), (0, 2), (0, 0)],
              [(0, 9), (1, 1), (0, 0), (2, 2), (1, 0), (4, 4)]]
    reports = [_report(tmp_path / f"r{i}.gz", c) for i, c in enumerate(counts)]
    fast, _ = _merge(tmp_path, reports)
    # Dropping one site from a report forces the sort path; add it back as a
    # separate one-line report so the totals stay the same.
    shifted = [_report(tmp_path / "s0.gz", counts[0][1:], SITES[1:]),
               _report(tmp_path / "s1.gz", counts[0][:1], SITES[:1]),
               *reports[1:]]
    slow, stdout = _merge(tmp_path, shifted)
    assert "merging by sort" in stdout
    assert _by_site(fast) == _by_site(slow)


def test_reports_covering_different_sites_merge_by_sort(tmp_path):
    a = _report(tmp_path / "a.gz", [(1, 0), (2, 3)], [SITES[0], SITES[2]])
    b = _report(tmp_path / "b.gz", [(4, 4), (1, 1)], [SITES[2], SITES[5]])
    rows, stdout = _merge(tmp_path, [a, b])
    assert "merging by sort" in stdout
    merged = _by_site(rows)
    assert merged[("chr2", 3)][3:5] == ["1", "0"]
    assert merged[("chr2", 40)][3:5] == ["6", "7"]
    assert merged[("chr1", 9)][3:5] == ["1", "1"]


def test_single_replicate_passes_through(tmp_path):
    counts = [(1, 0), (0, 0), (2, 3), (0, 4), (5, 5), (1, 1)]
    rows, _ = _merge(tmp_path, [_report(tmp_path / "a.gz", counts)])
    assert rows == [[c, str(p), s, str(m), str(u), ctx, tri]
                    for (c, p, s, ctx, tri), (m, u) in zip(SITES, counts)]


def test_corrupt_report_fails(tmp_path):
    a = _report(tmp_path / "a.gz", [(1, 0)] * len(SITES))
    bad = tmp_path / "bad.gz"
    bad.write_bytes(Path(a).read_bytes()[:-6])  # truncated gzip trailer
    out = tmp_path / "merged.gz"
    result = subprocess.run(["bash", str(SCRIPT), "1", str(out), a, str(bad)],
                            capture_output=True, text=True,
                            env={**os.environ, "TMPDIR": str(tmp_path)})
    assert result.returncode != 0
