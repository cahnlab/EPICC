"""matrix_io.R: reading computeMatrix files in R for the R plotting scripts.

The R side is checked against small matrices written here in deeptools'
format (a stranded group, a 'nostrand' group, NaNs, an all-NaN row, tied
rows, an empty group), with expected values computed in Python the way
deeptools computes them (np.nanmean / np.nanmedian, numpy argsort reversed
for 'descend'). Where deeptools' plotProfile and plotHeatmap are on PATH,
group_means and sort_rows are also compared with what deeptools itself writes
for the same matrix.
"""

import gzip
import json
import math
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
MATRIX_IO = REPO / "workflow" / "scripts" / "matrix_io.R"


def _have_r():
    if not shutil.which("Rscript"):
        return False
    proc = subprocess.run(["Rscript", "-e", "library(data.table); library(jsonlite)"],
                          capture_output=True)
    return proc.returncode == 0


needs_r = pytest.mark.skipif(not _have_r(), reason="Rscript with data.table/jsonlite not found")
needs_deeptools = pytest.mark.skipif(
    not (shutil.which("plotProfile") and shutil.which("plotHeatmap")),
    reason="deeptools not found")

NAN = float("nan")


def write_matrix(path, groups, samples, rows, *, upstream=1000, downstream=1000,
                 body=1500, bin_size=100, ref_point=None, unscaled5=0, unscaled3=0):
    """groups: [(label, n)], samples: [label], rows: [(bed6, values)]."""
    nbins = len(rows[0][1]) // len(samples) if rows else 0
    gb = [0]
    for _, n in groups:
        gb.append(gb[-1] + n)
    per = lambda v: [v] * len(samples)  # noqa: E731
    header = {
        "upstream": per(upstream), "downstream": per(downstream), "body": per(body),
        "bin size": per(bin_size), "ref point": per(ref_point),
        "verbose": False, "bin avg type": "mean", "missing data as zero": False,
        "min threshold": None, "max threshold": None, "scale": 1, "skip zeros": False,
        "nan after end": False, "proc number": 1, "sort regions": "keep",
        "sort using": "mean", "unscaled 5 prime": per(unscaled5),
        "unscaled 3 prime": per(unscaled3),
        "group_labels": [g for g, _ in groups], "group_boundaries": gb,
        "sample_labels": samples,
        "sample_boundaries": [i * nbins for i in range(len(samples) + 1)],
    }
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header, separators=(",", ":")) + "\n")
        for bed, vals in rows:
            fh.write("\t".join(map(str, bed)) + "\t"
                     + "\t".join("nan" if math.isnan(v) else f"{v:.6f}" for v in vals) + "\n")
    return header


def fixture_rows(nbins=6):
    """A stranded group (7 rows) and a nostrand group (4 rows), 2 samples.

    r2 is all NaN, r4 and r5 are identical (a tie), q1 has NaN in sample A,
    r6's sample B is all NaN."""
    rows = []
    for i in range(7):
        a = [((i * 7 + k * 3) % 11) / 2 for k in range(nbins)]
        b = [((i * 5 + k * 2) % 9) / 3 for k in range(nbins)]
        if i == 2:
            a = b = [NAN] * nbins
        if i == 6:
            b = [NAN] * nbins
        rows.append((("chr1", 100 * i, 100 * i + 50, f"r{i}", 0, "+-"[i % 2]), a + b))
    rows[5] = (rows[5][0], list(rows[4][1]))
    for i in range(4):
        a = [(i + k) / 4 for k in range(nbins)]
        b = [(i * k % 5) / 2 for k in range(nbins)]
        if i == 1:
            a = [NAN] * 3 + a[3:]
        rows.append((("chr2", 100 * i, 100 * i + 50, f"q{i}", 0, "."), a + b))
    return rows


def nan_summary(vals, how):
    v = [x for x in vals if not math.isnan(x)]
    if not v:
        return 0.0 if how == "sum" else NAN
    if how == "mean":
        return sum(v) / len(v)
    if how == "median":
        v = sorted(v)
        n = len(v)
        return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2
    return {"max": max, "min": min, "sum": sum}[how](v)


def expected_order(rows, groups, how, method):
    """deeptools' order: ascending (NaN last, stable), reversed for descend."""
    out, start = [], 0
    for _, n in groups:
        idx = list(range(start, start + n))
        score = [nan_summary(rows[i][1], how) for i in idx]
        asc = sorted(range(n), key=lambda j: (math.isnan(score[j]), 0 if math.isnan(score[j]) else score[j]))
        if method == "descend":
            asc = asc[::-1]
        out += [rows[idx[j]][0][3] for j in asc]
        start += n
    return out


def run_r(code):
    prog = f'source("{MATRIX_IO}")\n{code}'
    proc = subprocess.run(["Rscript", "-e", prog], capture_output=True, text=True)
    return proc


def run_r_ok(code):
    proc = run_r(code)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


GROUPS = [("all_TEs", 7), ("all_TEs_nostrand", 4)]


@pytest.fixture
def matrix(tmp_path):
    rows = fixture_rows()
    path = tmp_path / "m.gz"
    write_matrix(path, GROUPS, ["A", "B"], rows)
    return path, rows


@needs_r
class TestRead:
    def test_layout(self, matrix):
        path, rows = matrix
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            cat(dim(m$values), "\\n")
            cat(m$groups$label, "|", m$groups$first, "|", m$groups$last, "|", m$groups$n, "\\n")
            cat(m$samples$label, "|", m$samples$first, "|", m$samples$last, "|", m$samples$nbins, "\\n")
            cat(levels(m$regions$group), "|", as.character(m$regions$group[c(1, 8)]), "\\n")
            cat(m$regions$name, "\\n")
            cat(m$regions$strand[c(1, 2, 8)], class(m$regions$start), "\\n")
            cat(sum(is.na(m$values)), "\\n")
            cat(m$samples$upstream[1], m$samples$body[1], m$samples$bin_size[1], is.na(m$samples$ref_point[1]), "\\n")
        ''').splitlines()
        assert out[0].split() == ["11", "12"]
        assert out[1].split() == ["all_TEs", "all_TEs_nostrand", "|", "1", "8", "|", "7",
                                  "11", "|", "7", "4"]
        assert out[2].split() == ["A", "B", "|", "1", "7", "|", "6", "12", "|", "6", "6"]
        assert out[3].split() == ["all_TEs", "all_TEs_nostrand", "|", "all_TEs", "all_TEs_nostrand"]
        assert out[4].split() == [r[0][3] for r in rows]
        assert out[5].split() == ["+", "-", ".", "integer"]
        n_nan = sum(math.isnan(v) for _, vals in rows for v in vals)
        assert int(out[6]) == n_nan
        assert out[7].split() == ["1000", "1500", "100", "TRUE"]

    def test_values_match_file(self, matrix):
        path, rows = matrix
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            cat(format(m$values[1, ], digits = 10), "\\n")
            cat(sample_columns(m, "B"), "\\n")
        ''').splitlines()
        assert [float(x) for x in out[0].split()] == pytest.approx(rows[0][1])
        assert out[1].split() == [str(i) for i in range(7, 13)]

    def test_row_count_mismatch_is_an_error(self, tmp_path):
        rows = fixture_rows()
        path = tmp_path / "bad.gz"
        write_matrix(path, [("all_TEs", 12)], ["A", "B"], rows)
        proc = run_r(f'read_matrix("{path}")')
        assert proc.returncode != 0
        assert "announces 12 regions, file has 11" in proc.stderr

    def test_not_a_matrix(self, tmp_path):
        path = tmp_path / "x.gz"
        with gzip.open(path, "wt") as fh:
            fh.write("chr1\t0\t1\n")
        proc = run_r(f'read_matrix("{path}")')
        assert proc.returncode != 0
        assert "not a deeptools matrix" in proc.stderr

    def test_empty_group(self, tmp_path):
        rows = fixture_rows()[:7]
        path = tmp_path / "e.gz"
        write_matrix(path, [("all_TEs", 7), ("all_TEs_nostrand", 0)], ["A", "B"], rows)
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            cat(m$groups$n, "\\n")
            g <- group_means(m)
            cat(all(is.na(g$value[g$group == "all_TEs_nostrand"])), nrow(g), "\\n")
            cat(nrow(sort_rows(m)$values), "\\n")
        ''').splitlines()
        assert out[0].split() == ["7", "0"]
        assert out[1].split() == ["TRUE", "24"]
        assert out[2].strip() == "7"


@needs_r
class TestGroupMeans:
    @pytest.mark.parametrize("average", ["mean", "median", "max"])
    def test_nan_aware(self, matrix, average):
        path, rows = matrix
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            g <- group_means(m, "{average}")
            write.table(g, sep = "\\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        ''').splitlines()
        got = {}
        for line in out:
            s, grp, x, v = line.split("\t")
            got[(s, grp, int(x))] = NAN if v == "NA" else float(v)
        start = 0
        for grp, n in GROUPS:
            for si, s in enumerate(["A", "B"]):
                for x in range(6):
                    col = [rows[r][1][si * 6 + x] for r in range(start, start + n)]
                    exp = nan_summary(col, average)
                    assert got[(s, grp, x)] == pytest.approx(exp, nan_ok=True), (s, grp, x)
            start += n
        assert len(got) == 2 * 2 * 6


@needs_r
class TestSortRows:
    @pytest.mark.parametrize("how", ["mean", "median", "max", "sum"])
    @pytest.mark.parametrize("method", ["descend", "ascend"])
    def test_order(self, matrix, how, method):
        path, rows = matrix
        out = run_r_ok(f'''
            m <- sort_rows(read_matrix("{path}"), "{method}", "{how}")
            cat(m$regions$name, "\\n")
            cat(as.character(m$regions$group), "\\n")
            cat(format(m$values[, 1], digits = 10), "\\n")
        ''').splitlines()
        names = out[0].split()
        assert names == expected_order(rows, GROUPS, how, method)
        # rows stay within their group, and values travel with their region
        assert out[1].split() == ["all_TEs"] * 7 + ["all_TEs_nostrand"] * 4
        by_name = {r[0][3]: r[1][0] for r in rows}
        got = [NAN if v == "NA" else float(v) for v in out[2].split()]
        assert got == pytest.approx([by_name[n] for n in names], nan_ok=True)

    def test_all_nan_row_goes_first_under_descend(self, matrix):
        path, _ = matrix
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            cat(sort_rows(m, "descend")$regions$name[1], sort_rows(m, "ascend")$regions$name[7], "\\n")
        ''')
        assert out.split() == ["r2", "r2"]

    def test_tie_reversed_under_descend(self, matrix):
        path, _ = matrix
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            d <- sort_rows(m, "descend")$regions$name; a <- sort_rows(m, "ascend")$regions$name
            cat(match("r4", d) > match("r5", d), match("r4", a) < match("r5", a), "\\n")
        ''')
        assert out.split() == ["TRUE", "TRUE"]

    def test_keep(self, matrix):
        path, rows = matrix
        out = run_r_ok(f'cat(sort_rows(read_matrix("{path}"), "keep")$regions$name)')
        assert out.split() == [r[0][3] for r in rows]

    def test_using_samples(self, matrix):
        path, rows = matrix
        out = run_r_ok(f'cat(sort_rows(read_matrix("{path}"), "descend", "mean", samples = "B")$regions$name)')
        only_b = [(bed, vals[6:]) for bed, vals in rows]
        assert out.split() == expected_order(only_b, GROUPS, "mean", "descend")


@needs_r
class TestBinRows:
    def test_shrinks_with_nan_aware_means(self):
        out = run_r_ok('''
            v <- matrix(c(1, 3, NA, NA, 5, 7, 9, NA, 2, 4), ncol = 1)
            b <- bin_rows(v, 3)
            cat(dim(b), "|", format(b[, 1]), "|", attr(b, "rows"), "\\n")
            s <- bin_rows(v, 20)
            cat(dim(s), "|", attr(s, "rows"), "\\n")
        ''').split("\n")
        # rows 1-4 -> mean(1,3) = 2; rows 5-7 -> mean(5,7,9) = 7; rows 8-10 -> mean(2,4) = 3
        first = out[0].split("|")
        assert first[0].split() == ["3", "1"]
        assert [float(x) for x in first[1].split()] == [2, 7, 3]
        assert first[2].split() == ["4", "3", "3"]
        second = out[1].split("|")
        assert second[0].split() == ["10", "1"]
        assert second[1].split() == ["1"] * 10

    def test_all_nan_bin_stays_na(self):
        out = run_r_ok('''
            v <- matrix(c(NA, NA, 1, 2), ncol = 1)
            cat(bin_rows(v, 2)[, 1])
        ''')
        assert out.split() == ["NA", "1.5"]


@needs_r
class TestTicks:
    def test_scale_regions(self, matrix):
        path, _ = matrix
        out = run_r_ok(f'''
            t <- matrix_ticks(read_matrix("{path}"), 1, start_label = "start", end_label = "end")
            cat(t$x, "|", t$label, "\\n")
        ''')
        # deeptools getProfileTicks: 0, b, b+m-w, b+m+a-w, in bins of w
        assert out.split() == ["0", "10", "24", "34", "|", "-1.0", "start", "end", "1.0Kb"]

    def test_reference_point(self, tmp_path):
        rows = [(("chr1", 0, 10, "r0", 0, "+"), [1.0] * 20)]
        path = tmp_path / "rp.gz"
        write_matrix(path, [("genes", 1)], ["A"], rows, body=0, ref_point="TSS",
                     upstream=2000, downstream=2000, bin_size=200)
        out = run_r_ok(f'''
            t <- matrix_ticks(read_matrix("{path}"), 1, ref_label = "TSS")
            cat(t$x, "|", t$label, "\\n")
        ''')
        assert out.split() == ["0", "9.5", "19", "|", "-2.0", "TSS", "2.0Kb"]

    def test_unscaled_flanks(self, tmp_path):
        rows = [(("chr1", 0, 10, "r0", 0, "+"), [1.0] * 45)]
        path = tmp_path / "u.gz"
        write_matrix(path, [("genes", 1)], ["A"], rows, unscaled5=500, unscaled3=500)
        out = run_r_ok(f'''
            t <- matrix_ticks(read_matrix("{path}"), 1)
            cat(t$x, "|", paste0("[", t$label, "]"), "\\n")
        ''')
        assert out.split() == ["0", "10", "15", "30", "34", "44", "|",
                               "[-1.0]", "[TSS]", "[]", "[]", "[TES]", "[1.0Kb]"]


@needs_r
@needs_deeptools
class TestAgainstDeeptools:
    """group_means and sort_rows give what deeptools writes for the same matrix."""

    @pytest.mark.parametrize("average", ["mean", "median"])
    def test_group_means_vs_plotProfile(self, matrix, tmp_path, average):
        path, _ = matrix
        tab = tmp_path / f"p_{average}.tab"
        subprocess.run(["plotProfile", "-m", str(path), "-o", str(tmp_path / "p.png"),
                        "--averageType", average, "--outFileNameData", str(tab)],
                       check=True, capture_output=True)
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            g <- group_means(m, "{average}")
            tab <- read.delim("{tab}", header = FALSE, skip = 2, stringsAsFactors = FALSE,
                              na.strings = c("nan", "NA"))
            d <- 0
            for (r in seq_len(nrow(tab))) {{
              v <- as.numeric(tab[r, -(1:2)]); v <- v[seq_len(6)]
              ref <- g$value[g$sample == tab[r, 1] & g$group == tab[r, 2]]
              stopifnot(length(ref) == 6)
              d <- max(d, abs(ifelse(is.na(ref) & is.na(v), 0, ref - v)))
            }}
            cat(nrow(tab), d, "\\n")
        ''')
        n, diff = out.split()
        assert n == "4"
        assert float(diff) < 1e-9

    @pytest.mark.parametrize("how", ["mean", "median"])
    @pytest.mark.parametrize("method", ["descend", "ascend"])
    def test_sort_vs_plotHeatmap(self, matrix, tmp_path, how, method):
        path, _ = matrix
        bed = tmp_path / "sorted.bed"
        subprocess.run(["plotHeatmap", "-m", str(path), "-o", str(tmp_path / "h.png"),
                        "--sortRegions", method, "--sortUsing", how,
                        "--outFileSortedRegions", str(bed)],
                       check=True, capture_output=True)
        want = [line.split("\t")[3] for line in bed.read_text().splitlines()[1:]]
        out = run_r_ok(f'cat(sort_rows(read_matrix("{path}"), "{method}", "{how}")$regions$name)')
        assert out.split() == want
