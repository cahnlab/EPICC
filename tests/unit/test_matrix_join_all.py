"""The `all` env's matrix is the join of the `most` and mC matrices.

It used to be one computeMatrix run over every bigwig with the non-mC options
(`--missingDataAsZero --skipZeros`), which read an uncovered cytosine as 0 %
methylation: on the hg38 chr21 fixture the mean mCG over genes came out at
38 % instead of 72 %. `joining_all_matrix` now builds it from the two matrices
that are each right for their own data (matrix_join_all.py), and the phase 3
combined heatmap reads it, so this is the one implementation of that join.
"""

import gzip
import json
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
PROFILE = REPO / "profiles" / "default" / "config.yaml"
sys.path.insert(0, str(SCRIPTS))

from matrix_join_all import join_files  # noqa: E402

TOTAL = 800  # upstream 200 + body 400 + downstream 200


def _write(path, samples, bin_size, groups, upstream=200, downstream=200):
    """groups: [(label, [(name, start, values)])], values spanning every sample.

    Written the way computeMatrix writes: an '@'-prefixed JSON header, then one
    BED6 + values row per region, values as '%f' and 'nan'.
    """
    n = len(samples)
    nbins = TOTAL // bin_size
    header = {
        "upstream": [upstream] * n, "downstream": [downstream] * n,
        "body": [TOTAL - upstream - downstream] * n, "bin size": [bin_size] * n,
        "ref point": [None] * n, "verbose": False, "bin avg type": "mean",
        "missing data as zero": True, "min threshold": None,
        "max threshold": None, "scale": 1, "skip zeros": True,
        "nan after end": False, "proc number": 1, "sort regions": "keep",
        "sort using": "mean", "unscaled 5 prime": [0] * n, "unscaled 3 prime": [0] * n,
        "group_labels": [g for g, _ in groups], "group_boundaries": [0],
        "sample_labels": samples,
        "sample_boundaries": [i * nbins for i in range(n + 1)],
    }
    lines = []
    for _, regions in groups:
        header["group_boundaries"].append(header["group_boundaries"][-1] + len(regions))
        for name, start, values in regions:
            assert len(values) == nbins * n
            lines.append("\t".join(["chr1", str(start), str(start + 400), name, "0", "+"]
                                   + ["nan" if v is None else f"{v:f}" for v in values]))
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header) + "\n")
        fh.write("".join(l + "\n" for l in lines))


def _read(path):
    """(header, {group: [(name, [floats])]}) of a matrix file."""
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline().lstrip("@"))
        rows = [l.rstrip("\n").split("\t") for l in fh]
    b = header["group_boundaries"]
    return header, {g: [(r[3], [float(v) for v in r[6:]]) for r in rows[b[i]:b[i + 1]]]
                    for i, g in enumerate(header["group_labels"])}


def _join(tmp_path, most=None, mc=None):
    out = tmp_path / "all.gz"
    join_files(str(most) if most else None, str(mc) if mc else None, str(out))
    assert not list(tmp_path.glob("*.tmp")), "temporary body file left behind"
    return _read(out)


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
    h, groups = _join(tmp_path, *pair)
    assert h["sample_labels"] == ["WT__H3K9me2", "WT__RNAseq", "WT__mCG"]
    assert h["bin size"] == [50, 50, 50]
    assert h["sample_boundaries"] == [0, 16, 32, 48]
    assert h["upstream"] == [200] * 3 and h["body"] == [400] * 3
    assert h["group_labels"] == ["genes"] and h["group_boundaries"] == [0, 4]
    assert all(len(v) == 48 for _, v in groups["genes"])


def test_header_flags_are_booleans(pair, tmp_path):
    """computeMatrix always writes booleans; each half's fill is in the values."""
    h, _ = _join(tmp_path, *pair)
    assert h["missing data as zero"] is False
    assert h["skip zeros"] is False


def test_mc_bins_are_repeated_exactly(pair, tmp_path):
    _, groups = _join(tmp_path, *pair)
    assert groups["genes"][0][1][32:] == [10.0 + i // 2 for i in range(16)]


def test_uncovered_cytosines_stay_missing(pair, tmp_path):
    """The bug: the old `all` matrix read these as 0 % methylation."""
    _, groups = _join(tmp_path, *pair)
    gene_b = dict(groups["genes"])["geneB"]
    assert gene_b[:32] == [3.0] * 32
    assert all(math.isnan(v) for v in gene_b[32:])


def test_rows_skipped_by_skipzeros_are_zero_and_empty_rows_dropped(pair, tmp_path):
    _, groups = _join(tmp_path, *pair)
    names = [n for n, _ in groups["genes"]]
    # mC order (every region, target-file order); geneD has no data anywhere.
    assert names == ["geneA", "geneA", "geneB", "geneC"]
    assert groups["genes"][3][1] == [0.0] * 48


def test_duplicate_regions_matched_in_order(pair, tmp_path):
    _, groups = _join(tmp_path, *pair)
    (_, first), (_, second) = groups["genes"][:2]
    assert first[0] == 1.0 and first[32] == 10.0
    assert second[0] == 2.0 and second[32] == 20.0


def test_most_only_rows_keep_their_target_file_position(tmp_path):
    """A row the mC matrix lacks goes after its predecessor in `most`, not to
    the end of its group, so `--sortRegions keep` shows target-file order."""
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 100, [("genes", [
        ("x0", 500, [5.0] * 8),            # before every mC row
        ("g1", 1000, [1.0] * 8),
        ("x1", 1500, [6.0] * 8),           # between g1 and g3
        ("x2", 1600, [7.0] * 8),
        ("g3", 3000, [3.0] * 8),
    ])])
    _write(mc, ["s__mCG"], 100, [("genes", [
        ("g1", 1000, [50.0] * 8),
        ("g2", 2000, [60.0] * 8),          # not in most: no non-mC signal
        ("g3", 3000, [70.0] * 8),
        ("g4", 4000, [80.0] * 8),
    ])])
    _, groups = _join(tmp_path, most, mc)
    assert [n for n, _ in groups["genes"]] == ["x0", "g1", "x1", "x2", "g2", "g3", "g4"]
    x1 = dict(groups["genes"])["x1"]
    assert x1[0] == 6.0 and all(math.isnan(v) for v in x1[8:])


def test_groups_join_by_label(tmp_path):
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 100, [
        ("genes", [("g1", 0, [1.0] * 8)]),
        ("genes_nostrand", [("n1", 5000, [5.0] * 8)]),
    ])
    _write(mc, ["s__mCG"], 100, [
        ("genes", [("g1", 0, [50.0] * 8)]),
        ("genes_nostrand", [("n1", 5000, [60.0] * 8), ("n2", 6000, [70.0] * 8)]),
    ])
    h, groups = _join(tmp_path, most, mc)
    assert h["group_labels"] == ["genes", "genes_nostrand"]
    assert h["group_boundaries"] == [0, 1, 3]
    assert dict(groups["genes_nostrand"])["n2"][:8] == [0.0] * 8


def test_common_divisor_when_neither_bin_size_divides_the_other(tmp_path):
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 40, [("genes", [("g", 0, [float(i) for i in range(20)])])])
    _write(mc, ["s__mCG"], 100, [("genes", [("g", 0, [float(i) for i in range(8)])])])
    h, groups = _join(tmp_path, most, mc)
    assert h["bin size"] == [20, 20]
    assert h["sample_boundaries"] == [0, 40, 80]
    values = groups["genes"][0][1]
    assert values[:40] == [float(i // 2) for i in range(40)]
    assert values[40:] == [float(i // 5) for i in range(40)]


def test_layout_mismatch_is_an_error(tmp_path):
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 100, [("genes", [("g", 0, [1.0] * 8)])])
    _write(mc, ["s__mCG"], 100, [("genes", [("g", 0, [1.0] * 8)])], upstream=300, downstream=100)
    with pytest.raises(ValueError, match="upstream"):
        join_files(str(most), str(mc), str(tmp_path / "all.gz"))


def test_out_of_order_most_is_still_joined(tmp_path):
    """Not the usual case (computeMatrix keeps target-file order in both), but
    the reader holds skipped rows rather than mismatching them."""
    most = tmp_path / "most.gz"
    mc = tmp_path / "mc.gz"
    _write(most, ["s__ChIP"], 100, [("genes", [("g2", 2000, [2.0] * 8), ("g1", 1000, [1.0] * 8)])])
    _write(mc, ["s__mCG"], 100, [("genes", [("g1", 1000, [10.0] * 8), ("g2", 2000, [20.0] * 8)])])
    _, groups = _join(tmp_path, most, mc)
    rows = dict(groups["genes"])
    assert rows["g1"][0] == 1.0 and rows["g1"][8] == 10.0
    assert rows["g2"][0] == 2.0 and rows["g2"][8] == 20.0


@pytest.mark.parametrize("which", ["most", "mc"])
def test_single_input_keeps_every_row(pair, tmp_path, which):
    """With one input nothing is dropped, so `all` shows the same regions as
    that env (the mC matrix keeps all-NaN rows; no --skipZeros there)."""
    most, mc = pair
    out = tmp_path / "all.gz"
    subprocess.run([sys.executable, str(SCRIPT), "-o", str(out), f"--{which}",
                    str(most if which == "most" else mc)], check=True, capture_output=True)
    h, groups = _read(out)
    names = [n for n, _ in groups["genes"]]
    if which == "most":
        assert h["sample_labels"] == ["WT__H3K9me2", "WT__RNAseq"]
        assert names == ["geneA", "geneA", "geneB"]
    else:
        assert h["bin size"] == [100]
        assert names == ["geneA", "geneA", "geneB", "geneC", "geneD"]


# --- the rule ----------------------------------------------------------------

def test_all_env_is_joined_not_computed():
    text = SMK.read_text()
    assert re.search(r"^ruleorder:\s*joining_all_matrix\s*>\s*merging_matrix", text, re.M)
    block = text[text.index("rule joining_all_matrix:"):]
    block = block[:block.index("\nrule ", 1)]
    assert "matrix_join_all.py" in block
    assert not re.search(r"^\s*computeMatrix\b", block, re.M)


def test_rule_has_resources_and_threads():
    text = PROFILE.read_text()
    assert re.search(r"^\s+joining_all_matrix:\s*\*\w+", text, re.M)
    assert re.search(r"^\s+joining_all_matrix:\s*1\s*$", text, re.M)


# --- against deeptools and matrix_io.R ---------------------------------------

def _deeptools():
    return all(shutil.which(t) for t in
               ("computeMatrix", "computeMatrixOperations", "plotHeatmap", "plotProfile"))


@pytest.mark.requires_deeptools
@pytest.mark.skipif(not _deeptools(), reason="deeptools not on PATH")
def test_against_computematrix(tmp_path):
    """A bigwig with an uncovered region: the old `all` options give 0, the
    join NaN; and every deeptools tool the pipeline points at the matrix
    reads it."""
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
    join_files(str(most), str(mc), str(new))

    _, old_rows = _read(old)
    _, new_rows = _read(new)
    old_uncovered = dict(old_rows["r"])["uncovered"]
    new_uncovered = dict(new_rows["r"])["uncovered"]
    assert old_uncovered[40:] == [0.0] * 40          # the bug
    assert all(math.isnan(v) for v in new_uncovered[40:])
    assert dict(new_rows["r"])["covered"][40:] == [80.0] * 40
    assert new_uncovered[:40] == old_uncovered[:40]  # non-mC unchanged

    sorted_bed = tmp_path / "sorted.bed"
    for cmd in (["plotHeatmap", "-m", new, "-out", tmp_path / "h.png",
                 "--outFileSortedRegions", sorted_bed],
                ["plotProfile", "-m", new, "-out", tmp_path / "p.png", "--perGroup"],
                ["computeMatrixOperations", "relabel", "-m", new, "--groupLabels", "x",
                 "-o", tmp_path / "relabelled.gz"],
                ["computeMatrixOperations", "sort", "-m", new,
                 "-R", sorted_bed, "-o", tmp_path / "sorted.gz"]):
        subprocess.run([str(c) for c in cmd], check=True, capture_output=True)


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not on PATH")
def test_matrix_io_reads_the_join(pair, tmp_path):
    """matrix_io.R (data.table fread on `gzip -dc`) reads the two gzip members
    as one stream, with NaN in the mC columns."""
    _join(tmp_path, *pair)
    code = (
        f'source("{SCRIPTS / "matrix_io.R"}"); m <- read_matrix("{tmp_path / "all.gz"}"); '
        'cat(nrow(m$values), ncol(m$values), sum(is.na(m$values)), "\\n")'
    )
    res = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True)
    if res.returncode and "there is no package called" in res.stderr:
        pytest.skip("data.table not installed")
    assert res.returncode == 0, res.stderr
    # 4 rows x 48 columns; geneB's 16 mC values are the only NaN.
    assert res.stdout.split() == ["4", "48", "16"]
