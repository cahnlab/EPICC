"""plot_scales.R: per-column heatmap and profile limits, ported from the bash of
computing_matrix_scales, so the R plots keep `heatmaps_scales: type|sample|
default` when deeptools stops drawing them.

The fixture matrix has the cases that broke before or are easy to break:
  - two marks with ranges two orders of magnitude apart (H3K9me2, H3K4me3);
  - H3K9me next to H3K9me2 (a mark must not take a longer mark's columns);
  - a stranded track whose plus/minus columns interleave by sample
    (WT__RNA_plus, WT__RNA_minus, ddm1__RNA_plus, ddm1__RNA_minus), the case
    of 51d0cbb;
  - mC columns with NaNs (dataRange gives NaN -> the heat_m* fallback), an
    all-zero mC context and an all-zero non-mC mark;
  - an all-NaN column (no profile -> 0..0.01).

Expected values are computed here the way the bash rule does it. Where
deeptools is on PATH, the bash rule itself (dev/scales_rule.py) is run on the
same matrix and compared with R, for every mode and both averages.
"""

import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from .test_matrix_io import needs_r, write_matrix

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"
sys.path.insert(0, str(REPO / "dev"))
from scales_rule import parse_flags, run_scales_rule  # noqa: E402

needs_deeptools = pytest.mark.skipif(
    not (shutil.which("plotProfile") and shutil.which("computeMatrixOperations")),
    reason="deeptools not found")

NAN = float("nan")
NBINS = 5
MC = (100, 50, 20)
KEYS = [
    "WT__H3K9me2", "ddm1__H3K9me2",
    "WT__H3K9me",
    "WT__H3K4me3",
    "WT__RNA_plus", "WT__RNA_minus", "ddm1__RNA_plus", "ddm1__RNA_minus",
    "WT__mCG", "WT__mCHH",
    "WT__H3K27me3",
    "WT__empty",
]
GROUPS = [("genes", 6), ("genes_nostrand", 3)]


def column_values(key, row):
    """NBINS values of one column for region `row`."""
    k = [((row * 7 + b * 3 + len(key)) % 10) / 10 for b in range(NBINS)]
    mark = key.split("__")[-1]
    if mark == "H3K9me2":
        base = 40 if key.startswith("WT") else 20
        return [base + 10 * x for x in k]
    if mark == "H3K9me":
        return [5 + x for x in k]
    if mark == "H3K4me3":
        return [0.02 * x for x in k]
    if mark.startswith("RNA"):
        scale = 8 if mark.endswith("plus") else 2
        return [scale * x * (2 if key.startswith("ddm1") else 1) for x in k]
    if mark == "mCG":
        v = [60 + 20 * x for x in k]
        if row % 3 == 0:
            v[1] = NAN
        return v
    if mark in ("mCHH", "H3K27me3"):
        return [0.0] * NBINS
    if mark == "empty":
        return [NAN] * NBINS
    raise AssertionError(key)


def build(path):
    nrows = sum(n for _, n in GROUPS)
    rows = []
    for r in range(nrows):
        vals = [v for key in KEYS for v in column_values(key, r)]
        rows.append((("chr1", r * 100, r * 100 + 50, f"g{r}", 0, "+-."[r % 3]), vals))
    write_matrix(path, GROUPS, KEYS, rows, body=1000, bin_size=200,
                 upstream=0, downstream=0)
    return rows


def percentile(v, q):
    v = sorted(v)
    h = (len(v) - 1) * q
    lo = math.floor(h)
    return v[lo] + (h - lo) * (v[min(lo + 1, len(v) - 1)] - v[lo])


def nan_mean(v):
    v = [x for x in v if not math.isnan(x)]
    return sum(v) / len(v) if v else NAN


def pad(vals):
    v = [x for x in vals if not math.isnan(x)]
    if not v:
        return (0.0, 0.01)
    lo, hi = min(v), max(v)
    p = (hi - lo) * 0.1
    if p <= 0:
        p = abs(hi) * 0.1
    if p <= 0:
        p = 0.01
    return (lo - p, hi + p)


def expected(rows, mode):
    """What computing_matrix_scales writes, column by column."""
    cols = {}
    for i, key in enumerate(KEYS):
        values = [v for _, vals in rows for v in vals[i * NBINS:(i + 1) * NBINS]]
        rng = (NAN, NAN) if any(math.isnan(v) for v in values) else \
            (percentile(values, 0.1), percentile(values, 0.9))
        prof, start = [], 0
        for _, n in GROUPS:
            for b in range(NBINS):
                prof.append(nan_mean([rows[r][1][i * NBINS + b] for r in range(start, start + n)]))
            start += n
        cols[key] = (rng, prof)
    fallback = {"mCG": MC[0], "mCHG": MC[1], "mCHH": MC[2]}
    units = {}
    for key in KEYS:
        units.setdefault(key.split("__")[-1] if mode == "type" else key, []).append(key)
    out = {}
    for unit, keys in units.items():
        lo = [cols[k][0][0] for k in keys if not math.isnan(cols[k][0][0])]
        hi = [cols[k][0][1] for k in keys if not math.isnan(cols[k][0][1])]
        z = (min(lo) if lo else 0.0, max(hi) if hi else 0.0)
        if z == (0.0, 0.0):
            z = (0.0, float(fallback.get(keys[0].split("__")[-1], 0.005)))
        y = pad([v for k in keys for v in cols[k][1]])
        for k in keys:
            out[k] = z + y
    return [out[k] for k in KEYS]


def r_scales(matrix, mode, average="mean"):
    prog = f'''
        source("{SCRIPTS / "matrix_io.R"}"); source("{SCRIPTS / "plot_scales.R"}")
        s <- matrix_scales(read_matrix("{matrix}"), "{mode}", "{average}",
                           mc_zmax = c(mCG = {MC[0]}, mCHG = {MC[1]}, mCHH = {MC[2]}))
        write.table(s, sep = "\\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    '''
    proc = subprocess.run(["Rscript", "-e", prog], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    rows = [line.split("\t") for line in proc.stdout.splitlines()]
    return rows


@pytest.fixture
def matrix(tmp_path):
    path = tmp_path / "most.gz"
    return path, build(path)


@needs_r
class TestModes:
    @pytest.mark.parametrize("mode", ["type", "sample"])
    def test_matches_the_bash_rule_logic(self, matrix, mode):
        path, rows = matrix
        got = r_scales(path, mode)
        assert [r[0] for r in got] == KEYS                    # column order, one row each
        want = expected(rows, mode)
        for key, r, w in zip(KEYS, got, want):
            assert [float(x) for x in r[2:6]] == pytest.approx(w, rel=1e-9, abs=1e-12), key

    def test_type_mode_one_range_per_mark(self, matrix):
        path, _ = matrix
        got = {r[0]: [float(x) for x in r[2:6]] for r in r_scales(path, "type")}
        assert got["WT__H3K9me2"] == got["ddm1__H3K9me2"]
        # H3K9me is its own mark, not part of H3K9me2
        assert got["WT__H3K9me"] != got["WT__H3K9me2"]
        assert got["WT__H3K9me"][1] < 7 < got["WT__H3K9me2"][1]
        # marks two orders of magnitude apart keep their own ranges
        assert got["WT__H3K4me3"][1] < 0.05 and got["WT__H3K9me2"][1] > 40

    def test_interleaved_stranded_columns(self, matrix):
        path, _ = matrix
        got = r_scales(path, "type")
        by_key = {r[0]: [float(x) for x in r[2:6]] for r in got}
        # plus and minus are distinct marks; each column carries its own
        assert by_key["WT__RNA_plus"] == by_key["ddm1__RNA_plus"]
        assert by_key["WT__RNA_minus"] == by_key["ddm1__RNA_minus"]
        assert by_key["WT__RNA_plus"] != by_key["WT__RNA_minus"]
        marks = [r[1] for r in got]
        assert marks[4:8] == ["RNA_plus", "RNA_minus", "RNA_plus", "RNA_minus"]

    @pytest.mark.parametrize("mode", ["type", "sample"])
    def test_fallbacks(self, matrix, mode):
        path, _ = matrix
        got = {r[0]: [float(x) for x in r[2:6]] for r in r_scales(path, mode)}
        assert got["WT__mCG"][:2] == [0, MC[0]]      # NaNs -> dataRange NaN -> heat_mcg
        assert got["WT__mCHH"][:2] == [0, MC[2]]     # all-zero context -> heat_mchh
        assert got["WT__H3K27me3"][:2] == [0, 0.005]  # all-zero other mark
        assert got["WT__H3K27me3"][2:] == pytest.approx([-0.01, 0.01])  # flat zero profile
        assert got["WT__empty"][2:] == [0, 0.01]     # no profile at all

    def test_default_has_no_per_column_limits(self, matrix):
        path, _ = matrix
        got = r_scales(path, "default")
        assert len(got) == len(KEYS)
        assert all(r[2:6] == ["NA"] * 4 for r in got)


@needs_r
class TestLookups:
    def _run(self, matrix, code):
        prog = f'''
            source("{SCRIPTS / "matrix_io.R"}"); source("{SCRIPTS / "plot_scales.R"}")
            m <- read_matrix("{matrix}"); s <- matrix_scales(m, "type")
            {code}
        '''
        proc = subprocess.run(["Rscript", "-e", prog], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        return proc.stdout

    def test_scale_for_column_by_key_and_index(self, matrix):
        path, _ = matrix
        out = self._run(path, '''
            a <- scale_for_column(s, "ddm1__RNA_minus"); b <- scale_for_column(s, 8)
            c <- scale_for_column(s, "WT__RNA_plus")
            cat(identical(a, b), a$zmax == s$zmax[6], c$zmax == s$zmax[7], a$zmax != c$zmax, "\\n")
        ''')
        assert out.split() == ["TRUE", "TRUE", "TRUE", "TRUE"]

    def test_unknown_column_is_an_error(self, matrix):
        path, _ = matrix
        prog = f'''
            source("{SCRIPTS / "matrix_io.R"}"); source("{SCRIPTS / "plot_scales.R"}")
            s <- matrix_scales(read_matrix("{path}"), "type"); scale_for_column(s, "nope")
        '''
        proc = subprocess.run(["Rscript", "-e", prog], capture_output=True, text=True)
        assert proc.returncode != 0 and "No scale for column nope" in proc.stderr

    def test_shared_ylim_collapses_min_and_max(self, matrix):
        path, _ = matrix
        out = self._run(path, 'y <- shared_ylim(s); cat(y[1] == min(s$ymin), y[2] == max(s$ymax))')
        assert out.split() == ["TRUE", "TRUE"]

    def test_flags_in_column_order(self, matrix):
        path, _ = matrix
        out = self._run(path, 'f <- scale_flags(s); cat(f[["heatmap"]], "\\n", f[["profile"]], "\\n")')
        heat = parse_flags(out.splitlines()[0])
        assert len(heat["--zMin"]) == len(heat["--zMax"]) == len(KEYS)
        assert heat["--zMax"][8:10] == [MC[0], MC[2]]

    def test_default_scales(self, matrix):
        path, _ = matrix
        out = self._run(path, '''
            d <- default_scales(m)
            v <- m$values
            cat(d$zmin == min(v, na.rm = TRUE), d$zmax == max(v, na.rm = TRUE), d$ymin < d$ymax, "\\n")
        ''')
        # the fixture has NaNs, so deeptools' percentiles are NaN and the
        # colour range falls back to the data's own min/max
        assert out.split() == ["TRUE", "TRUE", "TRUE"]


@needs_r
@needs_deeptools
class TestAgainstTheRule:
    """plot_scales.R against the bash rule itself, on the same matrix."""

    @pytest.mark.parametrize("mode", ["type", "sample"])
    @pytest.mark.parametrize("average", ["mean", "median"])
    def test_same_flags(self, matrix, tmp_path, mode, average):
        path, _ = matrix
        heat, prof = run_scales_rule(path, tmp_path / "rule", mode, average, MC)
        got = r_scales(path, mode, average)
        want = {**parse_flags(heat), **parse_flags(prof)}
        for flag, col in (("--zMin", 2), ("--zMax", 3), ("--yMin", 4), ("--yMax", 5)):
            assert [float(r[col]) for r in got] == pytest.approx(want[flag], rel=1e-5, abs=1e-6), flag

    def test_default_writes_nothing(self, matrix, tmp_path):
        path, _ = matrix
        assert run_scales_rule(path, tmp_path / "rule", "default") == ("", "")
