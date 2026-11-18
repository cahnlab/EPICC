"""The `all` env's matrix is the join of the `most` and mC matrices.

It used to be one computeMatrix run over every bigwig with the non-mC options
(`--missingDataAsZero --skipZeros`, bin size `bs`), which read an uncovered
cytosine as 0 % methylation: on the hg38 chr21 fixture the mean mCG over genes
came out at 38 % instead of 72 %. `joining_all_matrix` now builds it from the
two matrices that are each right for their own data (matrix_join_all.py).
"""

import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"
SCRIPT = SCRIPTS / "matrix_join_all.py"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"
sys.path.insert(0, str(SCRIPTS))

from matrix_join_all import join, read_matrix, write_matrix  # noqa: E402


def _header(samples, bin_size, nbins, groups, **layout):
    n = len(samples)
    h = {
        "upstream": [layout.get("upstream", 200)] * n,
        "downstream": [layout.get("downstream", 200)] * n,
        "body": [layout.get("body", 400)] * n,
        "bin size": [bin_size] * n,
        "ref point": [None] * n,
        "verbose": False, "bin avg type": "mean",
        "missing data as zero": True, "min threshold": None,
        "max threshold": None, "scale": 1, "skip zeros": True,
        "nan after end": False, "proc number": 1, "sort regions": "keep",
        "sort using": "mean",
        "unscaled 5 prime": [0] * n, "unscaled 3 prime": [0] * n,
        "group_labels": [g for g, _ in groups],
        "group_boundaries": [0],
        "sample_labels": samples,
        "sample_boundaries": [i * nbins for i in range(n + 1)],
    }
    for _, rows in groups:
        h["group_boundaries"].append(h["group_boundaries"][-1] + len(rows))
    return h


def _write(path, samples, bin_size, groups, **layout):
    """groups: [(label, [(name, start, values)])], values spanning every sample."""
    total = 800  # upstream + body + downstream
    nbins = total // bin_size
    header = _header(samples, bin_size, nbins, groups, **layout)
    rows = []
    for _, regions in groups:
        for name, start, values in regions:
            assert len(values) == nbins * len(samples)
            rows.append(["chr1", str(start), str(start + 400), name, "0", "+"]
                        + [("nan" if v is None else f"{v:f}") for v in values])
    write_matrix(path, header, rows)


def _rows(path):
    header, groups = read_matrix(path)
    return header, {g: [(bed[3], [float(v) for v in vals]) for _, bed, vals in rows]
                    for g, rows in groups.items()}


@pytest.fixture
def pair(tmp_path):
    """most: 2 samples at bin 50 (16 bins each); mC: 1 sample at bin 100 (8 bins).

    most lacks geneC (--skipZeros dropped it: no non-mC signal). mC has every
    region, geneB uncovered (all NaN) and geneD uncovered too -- geneD is also
    absent from most, so it carries no data anywhere. geneA appears twice
    (duplicate name and coordinates), matched in order.
    """
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["WT__H3K9me2", "WT__RNAseq"], 50, [("genes", [
        ("geneA", 0, [1.0] * 32),
        ("geneA", 0, [2.0] * 32),
        ("geneB", 1000, [3.0] * 32),
    ])])
    _write(mc, ["WT__mCG"], 100, [("genes", [
        ("geneA", 0, [10.0 + i for i in range(8)]),
        ("geneA", 0, [20.0] * 8),
        ("geneB", 1000, [None] * 8),
        ("geneC", 2000, [0.0] * 8),
        ("geneD", 3000, [None] * 8),
    ])])
    return most, mc


def test_columns_header_and_common_bin_size(pair, tmp_path):
    most, mc = pair
    out = tmp_path / "all.gz"
    header, rows = join(read_matrix(most), read_matrix(mc))
    write_matrix(out, header, rows)
    h, groups = _rows(out)
    assert h["sample_labels"] == ["WT__H3K9me2", "WT__RNAseq", "WT__mCG"]
    assert h["bin size"] == [50, 50, 50]
    assert h["sample_boundaries"] == [0, 16, 32, 48]
    assert h["upstream"] == [200] * 3 and h["body"] == [400] * 3
    assert h["group_labels"] == ["genes"] and h["group_boundaries"] == [0, 4]
    assert all(len(v) == 48 for _, v in groups["genes"])


def test_mc_bins_are_repeated_exactly(pair):
    most, mc = pair
    _, rows = join(read_matrix(most), read_matrix(mc))
    mc_part = [float(v) for v in rows[0][6 + 32:]]
    assert mc_part == [10.0 + i // 2 for i in range(16)]


def test_uncovered_cytosines_stay_missing(pair):
    """The bug: the old `all` matrix read these as 0 % methylation."""
    most, mc = pair
    _, rows = join(read_matrix(most), read_matrix(mc))
    gene_b = [r for r in rows if r[3] == "geneB"][0]
    assert [float(v) for v in gene_b[6:38]] == [3.0] * 32
    assert all(math.isnan(float(v)) for v in gene_b[38:])


def test_rows_skipped_by_skipzeros_are_zero_and_empty_rows_dropped(pair):
    most, mc = pair
    _, rows = join(read_matrix(most), read_matrix(mc))
    names = [r[3] for r in rows]
    # mC order (every region, target-file order); geneD has no data anywhere.
    assert names == ["geneA", "geneA", "geneB", "geneC"]
    gene_c = rows[3]
    assert [float(v) for v in gene_c[6:]] == [0.0] * 48


def test_duplicate_regions_matched_in_order(pair):
    most, mc = pair
    _, rows = join(read_matrix(most), read_matrix(mc))
    first, second = rows[0], rows[1]
    assert float(first[6]) == 1.0 and float(first[38]) == 10.0
    assert float(second[6]) == 2.0 and float(second[38]) == 20.0


def test_groups_join_by_label_and_most_only_rows_are_kept(tmp_path):
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 100, [
        ("genes", [("g1", 0, [1.0] * 8), ("g9", 9000, [9.0] * 8)]),
        ("genes_nostrand", [("n1", 5000, [5.0] * 8)]),
    ])
    _write(mc, ["s__mCG"], 100, [
        ("genes", [("g1", 0, [50.0] * 8)]),
        ("genes_nostrand", [("n1", 5000, [60.0] * 8), ("n2", 6000, [70.0] * 8)]),
    ])
    header, rows = join(read_matrix(most), read_matrix(mc))
    assert header["group_labels"] == ["genes", "genes_nostrand"]
    assert header["group_boundaries"] == [0, 2, 4]
    assert [r[3] for r in rows] == ["g1", "g9", "n1", "n2"]
    g9 = rows[1]
    assert float(g9[6]) == 9.0 and all(math.isnan(float(v)) for v in g9[14:])
    n2 = rows[3]
    assert float(n2[6]) == 0.0 and float(n2[14]) == 70.0


def test_common_divisor_when_neither_bin_size_divides_the_other(tmp_path):
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 40, [("genes", [("g", 0, [float(i) for i in range(20)])])])
    _write(mc, ["s__mCG"], 100, [("genes", [("g", 0, [float(i) for i in range(8)])])])
    header, rows = join(read_matrix(most), read_matrix(mc))
    assert header["bin size"] == [20, 20]
    assert header["sample_boundaries"] == [0, 40, 80]
    values = [float(v) for v in rows[0][6:]]
    assert values[:40] == [float(i // 2) for i in range(40)]
    assert values[40:] == [float(i // 5) for i in range(40)]


def test_layout_mismatch_is_an_error(tmp_path):
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 100, [("genes", [("g", 0, [1.0] * 8)])])
    _write(mc, ["s__mCG"], 100, [("genes", [("g", 0, [1.0] * 8)])], upstream=300, downstream=100)
    with pytest.raises(ValueError, match="upstream"):
        join(read_matrix(most), read_matrix(mc))


@pytest.mark.parametrize("which", ["most", "mc"])
def test_single_input(pair, tmp_path, which):
    most, mc = pair
    out = tmp_path / "all.gz"
    subprocess.run([sys.executable, str(SCRIPT), "-o", str(out), f"--{which}",
                    str(most if which == "most" else mc)], check=True)
    h, groups = _rows(out)
    names = [n for n, _ in groups["genes"]]
    if which == "most":
        assert h["sample_labels"] == ["WT__H3K9me2", "WT__RNAseq"]
        assert names == ["geneA", "geneA", "geneB"]
    else:
        assert h["bin size"] == [100]
        # all-NaN rows carry no data
        assert names == ["geneA", "geneA", "geneC"]


# --- the rule ----------------------------------------------------------------

def test_all_env_is_joined_not_computed():
    text = SMK.read_text()
    assert re.search(r"^ruleorder:\s*joining_all_matrix\s*>\s*merging_matrix", text, re.M)
    block = text[text.index("rule joining_all_matrix:"):]
    block = block[:block.index("\nrule ", 1)]
    assert "matrix_join_all.py" in block
    assert not re.search(r"^\s*computeMatrix\b", block, re.M)


# --- against deeptools ------------------------------------------------------

def _deeptools():
    return all(shutil.which(t) for t in ("computeMatrix", "plotHeatmap", "plotProfile"))


@pytest.mark.requires_deeptools
@pytest.mark.skipif(not _deeptools(), reason="deeptools not on PATH")
def test_against_computematrix(tmp_path):
    """A bigwig with an uncovered region: the old `all` options give 0, the join NaN."""
    pyBigWig = pytest.importorskip("pyBigWig")
    for name, (start, end, value) in {"mc.bw": (0, 3000, 80.0), "chip.bw": (0, 10000, 2.0)}.items():
        bw = pyBigWig.open(str(tmp_path / name), "w")
        bw.addHeader([("chr1", 10000)])
        bw.addEntries(["chr1"], [start], ends=[end], values=[value])
        bw.close()
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t1000\t2000\tcovered\t0\t+\nchr1\t6000\t7000\tuncovered\t0\t+\n")
    common = ["-R", str(bed), "-b", "500", "-a", "500", "-m", "1000"]

    def cm(out, opts, bws, labels, bs):
        subprocess.run(["computeMatrix", "scale-regions", *opts, *common, "-S", *bws,
                        "--samplesLabel", *labels, "-bs", str(bs), "-o", str(out)],
                       check=True, capture_output=True)

    old = tmp_path / "old_all.gz"
    cm(old, ["--missingDataAsZero", "--skipZeros"], [tmp_path / "chip.bw", tmp_path / "mc.bw"],
       ["s__ChIP", "s__mCG"], 50)
    most = tmp_path / "most.gz"
    cm(most, ["--missingDataAsZero", "--skipZeros"], [tmp_path / "chip.bw"], ["s__ChIP"], 50)
    mc = tmp_path / "mc.gz"
    cm(mc, [], [tmp_path / "mc.bw"], ["s__mCG"], 100)
    new = tmp_path / "all.gz"
    subprocess.run([sys.executable, str(SCRIPT), "-o", str(new), "--most", str(most),
                    "--mc", str(mc)], check=True)

    _, old_rows = _rows(old)
    _, new_rows = _rows(new)
    old_uncovered = dict(old_rows["r"])["uncovered"]
    new_uncovered = dict(new_rows["r"])["uncovered"]
    assert old_uncovered[40:] == [0.0] * 40          # the bug
    assert all(math.isnan(v) for v in new_uncovered[40:])
    assert dict(new_rows["r"])["covered"][40:] == [80.0] * 40
    assert new_uncovered[:40] == old_uncovered[:40]  # non-mC unchanged

    for cmd in (["plotHeatmap", "-m", str(new), "-out", str(tmp_path / "h.png")],
                ["plotProfile", "-m", str(new), "-out", str(tmp_path / "p.png"), "--perGroup"]):
        subprocess.run(cmd, check=True, capture_output=True)
