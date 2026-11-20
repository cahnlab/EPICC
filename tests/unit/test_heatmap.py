"""R_heatmap.R / plot_heatmap.R: Heatmap__ drawn in R with ComplexHeatmap.

Checked without looking at pixels where possible: the resampler on small
matrices with known pixel values (both interpolation modes, NaN edges, up- and
downsampling), the row order against deeptools' rules (numpy argsort reversed
under descend, all-NaN rows on top; in the all env on the non-mC columns only),
z limits looked up per column by key (interleaved stranded columns), the
region groups (stranded + nostrand), and the script end to end (PDF and
--outFileSortedRegions-style BED). Where deeptools' plotHeatmap is on PATH,
the sorted BED is compared with the one it writes.
"""

import math
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from .test_matrix_io import needs_r, write_matrix

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"
SCRIPT = SCRIPTS / "R_heatmap.R"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"
NAN = float("nan")
NBINS = 4


def _have_complexheatmap():
    if not shutil.which("Rscript"):
        return False
    proc = subprocess.run(["Rscript", "-e", "library(ComplexHeatmap); library(data.table)"],
                          capture_output=True)
    return proc.returncode == 0


needs_ch = pytest.mark.skipif(not _have_complexheatmap(),
                              reason="Rscript with ComplexHeatmap not found")
needs_deeptools = pytest.mark.skipif(not shutil.which("plotHeatmap"), reason="deeptools not found")

# An all-env matrix: ChIP, an unstranded matrix's interleaved RNA plus/minus
# columns, then mC contexts; a stranded group and a nostrand group.
ALL_KEYS = ["WT__H3K9me2", "ddm1__H3K9me2", "WT__RNAseq_plus", "WT__RNAseq_minus",
            "WT__mCG", "ddm1__mCG", "WT__mCHH"]
TYPES = {"H3K9me2": "ChIP", "RNAseq_plus": "RNA", "RNAseq_minus": "RNA",
         "mCG": "mCG", "mCHH": "mCHH"}
GROUPS = [("all_genes", 6), ("all_genes_nostrand", 3)]
COLOURS = {"WT": "#1f77b4", "ddm1": "#ff7f0e"}


def value(key, r, b):
    """Row r's bin b of column key: distinct row means, mC with gaps."""
    mark = key.split("__")[1]
    v = ((r * 7 + 3) % 11) + b * 0.1
    if mark.startswith("mC"):
        if (r + b) % 4 == 0:
            return NAN
        # mC ordered against the ChIP rows, so a sort on mC would differ
        return 100 - v * 5
    if mark == "RNAseq_minus":
        return 50 - v
    return v * (2 if key.startswith("ddm1") else 1)


def build(path, keys=ALL_KEYS, groups=GROUPS, nan_rows=(), nbins=NBINS):
    rows, r = [], 0
    for g, n in groups:
        for i in range(n):
            vals = []
            for k in keys:
                vals += [NAN if r in nan_rows else value(k, r, b) for b in range(nbins)]
            strand = "." if g.endswith("nostrand") else "+-"[i % 2]
            rows.append((("chr1", r * 100, r * 100 + 50, f"{g}_{i}", 0, strand), vals))
            r += 1
    write_matrix(path, groups, keys, rows, body=nbins * 100, bin_size=100,
                 upstream=0, downstream=0)
    return rows


def columns_tsv(path, keys=ALL_KEYS, types=True):
    lines = []
    for k in keys:
        sample, track = k.split("__")
        row = [k, f"{sample} {track}", COLOURS[sample], track, sample]
        if types:
            row.append(TYPES[track])
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")


def scales_txt(path, keys=ALL_KEYS, empty=False):
    if empty:
        path.write_text("")
        return
    zmin = [str(i) for i in range(len(keys))]
    zmax = [str(10 + i) for i in range(len(keys))]
    path.write_text(f"--zMin {' '.join(zmin)} --zMax {' '.join(zmax)}")


def expected_order(rows, groups, cols, nbins=NBINS, using="mean"):
    """deeptools: per group, numpy argsort of the row summary (NaN last),
    reversed for descend; summaries NaN-aware over the given columns."""
    out, start = [], 0
    for _, n in groups:
        scores = []
        for r in range(start, start + n):
            v = [rows[r][1][c * nbins + b] for c in cols for b in range(nbins)]
            v = [x for x in v if not math.isnan(x)]
            if not v:
                s = NAN
            elif using == "mean":
                s = sum(v) / len(v)
            else:
                v = sorted(v)
                k = len(v)
                s = v[k // 2] if k % 2 else (v[k // 2 - 1] + v[k // 2]) / 2
            scores.append(s)
        idx = sorted(range(n), key=lambda i: (math.isnan(scores[i]), scores[i] if not math.isnan(scores[i]) else 0))
        out += [start + i for i in reversed(idx)]
        start += n
    return out


def run_r(code):
    prog = "".join(f'source("{SCRIPTS / f}")\n' for f in
                   ("matrix_io.R", "plot_scales.R", "plot_theme.R", "plot_profile.R",
                    "plot_heatmap.R")) + code
    return subprocess.run(["Rscript", "-e", prog], capture_output=True, text=True)


def run_r_ok(code):
    proc = run_r(code)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


def run_script(tmp, matrix, cols, scales, env="all", mode="type", sort="mean", theme=None):
    pdf, bed = tmp / "Heatmap.pdf", tmp / "sorted.bed"
    args = ["Rscript", str(SCRIPT), str(matrix), str(cols), "regions", env, mode, "mean", sort,
            str(scales), str(pdf), str(bed)]
    if theme is not None:
        args.append(str(theme))
    return subprocess.run(args, capture_output=True, text=True), pdf, bed


def bed_names(bed):
    lines = bed.read_text().splitlines()
    assert lines[0].startswith("#chrom\tstart\tend\tname\tscore\tstrand")
    return [ln.split("\t") for ln in lines[1:]]


# --- resampling ---------------------------------------------------------------

def _matrix(rows):
    flat = ", ".join("NA" if (isinstance(x, float) and math.isnan(x)) else repr(float(x))
                     for row in rows for x in row)
    return f"matrix(c({flat}), nrow = {len(rows)}, byrow = TRUE)"


def _out(code):
    return [NAN if x == "NA" else float(x) for x in run_r_ok(code).split()]


@needs_ch
class TestResample:
    def test_nearest_repeats_rows(self):
        got = _out(f"cat(resample_rows({_matrix([[1], [2], [3]])}, 6, 'nearest'))")
        assert got == [1, 1, 2, 2, 3, 3]

    def test_bilinear_between_row_centres_edges_clamped(self):
        got = _out(f"cat(resample_rows({_matrix([[0], [10]])}, 4, 'bilinear'))")
        assert got == pytest.approx([0, 2.5, 7.5, 10])

    def test_bilinear_never_mixes_in_a_missing_value(self):
        got = _out(f"cat(resample_rows({_matrix([[0], [NAN], [10]])}, 6, 'bilinear'))")
        assert got[:2] == [0, 0] and got[4:] == [10, 10]
        assert all(math.isnan(x) for x in got[2:4])

    def test_downsampling_averages_runs_nan_aware(self):
        for mode in ("bilinear", "nearest"):
            got = _out(f"cat(resample_rows({_matrix([[1], [3], [NAN], [5]])}, 2, '{mode}'))")
            assert got == [2, 5]

    def test_same_size_is_unchanged(self):
        got = _out(f"cat(resample_rows({_matrix([[1, 2], [3, NAN]])}, 2, 'bilinear'))")
        assert got[:3] == [1, 3, 2] and math.isnan(got[3])

    def test_image_resamples_bins_too(self):
        # 1 row x 2 bins -> 2 x 4 pixels.
        bil = _out(f"cat(t(resample_image({_matrix([[0, 10]])}, 2, 4, 'bilinear')))")
        assert bil == pytest.approx([0, 2.5, 7.5, 10] * 2)
        near = _out(f"cat(t(resample_image({_matrix([[0, 10]])}, 2, 4, 'nearest')))")
        assert near == [0, 0, 10, 10] * 2

    def test_missing_bin_stays_a_block(self):
        # A missing mC bin between two covered ones, at 3x: 3 missing pixels.
        got = _out(f"cat(resample_image({_matrix([[20, NAN, 80]])}, 1, 9, 'bilinear'))")
        assert [math.isnan(x) for x in got] == [False] * 3 + [True] * 3 + [False] * 3
        assert got[:2] == [20, 20] and got[-2:] == [80, 80]

    def test_colours_clamp_and_mark_missing(self):
        out = run_r_ok("f <- heatmap_colour_fun(c('#000000', '#FFFFFF'), 0, 10, 'grey90'); "
                       "cat(f(c(-5, 0, 5, 10, 99, NA)))")
        assert out.split() == ["#000000", "#000000", "#808080", "#FFFFFF", "#FFFFFF", "#E5E5E5"]

    @pytest.mark.parametrize("z", [0, 5, -3])
    def test_a_flat_range_is_drawn_mid_ramp(self, z):
        # A column with one value throughout (an sRNA size class with no
        # reads) takes the ramp's middle colour, white on a diverging map, not
        # its low end, which would read as strong depletion.
        out = run_r_ok(f"f <- heatmap_colour_fun(c('#0000FF', '#FFFFFF', '#FF0000'), {z}, {z}, 'grey90'); "
                       f"cat(f({z}), range(attr(f, 'breaks'))[1] < {z})")
        assert out.split() == ["#FFFFFF", "TRUE"]


# --- order, scales, groups ----------------------------------------------------

@needs_ch
class TestPrepare:
    @pytest.fixture
    def all_matrix(self, tmp_path):
        path, cols = tmp_path / "all.gz", tmp_path / "cols.tsv"
        rows = build(path, nan_rows=(1, 7))
        columns_tsv(cols)
        return path, cols, rows

    @pytest.mark.parametrize("using", ["mean", "median"])
    def test_row_order_matches_deeptools_rules(self, all_matrix, using):
        path, _, rows = all_matrix
        out = run_r_ok(f'm <- read_matrix_table("{path}"); '
                       f'cat(row_order(m, "descend", "{using}") - 1)')
        got = [int(x) for x in out.split()]
        assert got == expected_order(rows, GROUPS, range(len(ALL_KEYS)), using=using)
        # all-NaN rows go to the top of their group under descend
        assert got[0] == 1 and got[6] == 7

    def test_table_and_matrix_readers_agree(self, all_matrix):
        path, _, _ = all_matrix
        out = run_r_ok(f'a <- read_matrix("{path}"); b <- read_matrix_table("{path}"); '
                       'cat(identical(row_order(a, "descend", "mean"), row_order(b, "descend", "mean")), '
                       'identical(matrix_block(a, 3), matrix_block(b, 3)), '
                       'identical(sort_rows(a)$values[, 1:4], matrix_block(sort_rows(b), 1)))')
        assert out.split() == ["TRUE", "TRUE", "TRUE"]

    def test_all_env_sorts_on_the_non_mc_columns(self, all_matrix):
        path, cols, rows = all_matrix
        out = run_r_ok(f'd <- read_display_table("{cols}"); cat(sort_columns(d, "all"), "|", '
                       f'is.null(sort_columns(d, "most")))')
        assert out.split() == ALL_KEYS[:4] + ["|", "TRUE"]
        out = run_r_ok(f'm <- read_matrix_table("{path}"); d <- read_display_table("{cols}"); '
                       f'cat(row_order(m, "descend", "mean", sort_columns(d, "all")) - 1)')
        assert [int(x) for x in out.split()] == expected_order(rows, GROUPS, range(4))
        # ... which is not the order a sort on every column gives here
        assert expected_order(rows, GROUPS, range(4)) != \
            expected_order(rows, GROUPS, range(len(ALL_KEYS)))

    def test_mc_only_sorts_on_mc(self, tmp_path):
        keys = ["WT__mCG", "ddm1__mCG"]
        cols = tmp_path / "c.tsv"
        columns_tsv(cols, keys)
        out = run_r_ok(f'cat(is.null(sort_columns(read_display_table("{cols}"), "all")))')
        assert out.strip() == "TRUE"

    def test_types_fall_back_on_the_mark_without_a_type_column(self, tmp_path):
        cols = tmp_path / "c.tsv"
        columns_tsv(cols, types=False)
        out = run_r_ok(f'cat(column_types(read_display_table("{cols}")))')
        assert out.split() == ["default"] * 4 + ["mCG", "mCG", "mCHH"]

    def test_z_limits_by_key(self, tmp_path):
        scales = tmp_path / "s.txt"
        scales_txt(scales)
        keys = ", ".join(f'"{k}"' for k in ALL_KEYS)
        out = run_r_ok(f'z <- read_heatmap_zlim("{scales}", c({keys})); '
                       'cat(unlist(scale_for_column(z, "WT__RNAseq_minus")[c("zmin", "zmax")]), '
                       'unlist(scale_for_column(z, "WT__mCHH")[c("zmin", "zmax")]))')
        assert [float(x) for x in out.split()] == [3, 13, 6, 16]
        scales.write_text("")
        assert run_r_ok(f'cat(is.null(read_heatmap_zlim("{scales}", c({keys}))))').strip() == "TRUE"
        scales.write_text("--zMin 0 1 --zMax 2 3")
        proc = run_r(f'read_heatmap_zlim("{scales}", c({keys}))')
        assert proc.returncode != 0 and "one per matrix column" in proc.stderr

    def test_default_scale_column_by_column_equals_one_vector(self):
        # heatmaps_scales: default takes the 1st/98th percentiles of every
        # value; R_heatmap.R works them out without gathering the values.
        out = run_r_ok("""
            set.seed(7)
            ok <- TRUE
            for (t in 1:60) {
              n <- sample(c(1, 2, 3, 10, 400), 1)
              cols <- lapply(1:sample(1:5, 1), function(j) switch(sample(3, 1),
                round(rnorm(n), 1), rep(3, n), sample(c(0, 0, 1, 9), n, TRUE)))
              if (t %% 10 == 0) cols[[1]][1] <- NA
              ok <- ok && identical(default_zlim(unlist(cols)),
                                    default_zlim_columns(function(j) cols[[j]], length(cols)))
            }
            cat(ok)""")
        assert out.strip() == "TRUE"

    def test_default_scale_per_data_type(self, all_matrix):
        # heatmaps_scales: default in the all env: each data type (ChIP, RNA,
        # mCG, mCHH here) gets a range of its own, not one shared.
        path, cols, rows = all_matrix
        out = run_r_ok(
            f'm <- read_matrix_table("{path}"); d <- read_display_table("{cols}"); '
            'z <- default_heatmap_zlim(m, d); '
            'v <- function(k) unlist(lapply(k, function(s) unlist(m$table[, sample_columns(m, s), with = FALSE]))); '
            'g <- list(1:2, 3:4, 5:6, 7); '
            'want <- do.call(rbind, lapply(g, function(k) t(replicate(length(k), default_zlim(v(k)))))); '
            'cat(identical(z$zmin, want[, 1]), identical(z$zmax, want[, 2]), '
            'z$zmax[1] != z$zmax[3] && z$zmax[1] < z$zmin[5])')
        assert out.split() == ["TRUE", "TRUE", "TRUE"]

    def test_default_scale_checks_the_display_first(self, all_matrix, tmp_path):
        # A display table that does not match the matrix gets the clear
        # message, not an indexing error from the per-type grouping.
        path, cols, rows = all_matrix
        short = tmp_path / "short.tsv"
        short.write_text("\n".join(cols.read_text().splitlines()[:-1]) + "\n")
        proc = run_r(f'm <- read_matrix_table("{path}"); '
                     f'default_heatmap_zlim(m, read_display_table("{short}"))')
        assert proc.returncode != 0
        assert "Display rows do not match the matrix columns" in proc.stderr

    def test_images_per_group_in_row_order(self, all_matrix):
        path, cols, rows = all_matrix
        out = run_r_ok(
            f'm <- read_matrix_table("{path}"); d <- read_display_table("{cols}"); th <- epicc_theme(); '
            'p <- prepare_heatmap(m, d, th, "descend", "mean", "mean", "most"); '
            'cat(length(p$images), length(p$images[[1]]), nrow(p$images[[1]][[1]]), '
            'nrow(p$images[[1]][[2]]), "|", p$images[[1]][[1]][, 1], "|", p$ord - 1)')
        head, first, order = out.split("|")
        assert head.split() == ["7", "2", "6", "3"]
        order = [int(x) for x in order.split()]
        assert order == expected_order(rows, GROUPS, range(len(ALL_KEYS)))
        want = [rows[r][1][0] for r in order[:6]]
        got = [NAN if x == "NA" else float(x) for x in first.split()]
        assert len(got) == 6
        for g, w in zip(got, want):
            assert (math.isnan(g) and math.isnan(w)) or g == pytest.approx(w)

    def test_many_regions_are_binned_to_the_pixel_rows(self, tmp_path):
        path, cols = tmp_path / "big.gz", tmp_path / "c.tsv"
        groups = [("genes", 400), ("genes_nostrand", 100)]
        build(path, keys=ALL_KEYS[:2], groups=groups)
        columns_tsv(cols, ALL_KEYS[:2])
        # 1 in x 100 dpi = 100 pixel rows for 500 regions: 80 + 20.
        out = run_r_ok(
            f'm <- read_matrix_table("{path}"); d <- read_display_table("{cols}"); '
            'th <- epicc_theme(); th$heatmap_height <- 1; th$heatmap_raster_dpi <- 100; '
            'p <- prepare_heatmap(m, d, th, "keep", "mean", "mean", "most"); '
            'cat(sapply(p$images[[2]], nrow))')
        assert out.split() == ["80", "20"]

    def test_profiles_are_profile_data_of_every_column(self, all_matrix):
        path, cols, _ = all_matrix
        out = run_r_ok(
            f'm <- read_matrix_table("{path}"); d <- read_display_table("{cols}"); th <- epicc_theme(); '
            'p <- prepare_heatmap(m, d, th, "descend", "mean", "mean", "all"); '
            f'full <- profile_data(read_matrix("{path}"), "mean"); '
            'cat(isTRUE(all.equal(p$d[, c("x", "value")], full[, c("x", "value")], check.attributes = FALSE)), '
            'identical(as.character(p$d$sample), as.character(full$sample)))')
        assert out.split() == ["TRUE", "TRUE"]


# --- the script ---------------------------------------------------------------

@needs_ch
class TestScript:
    @pytest.mark.parametrize("mode", ["type", "sample", "default"])
    def test_all_env_end_to_end(self, tmp_path, mode):
        path, cols, scales = tmp_path / "all.gz", tmp_path / "cols.tsv", tmp_path / "s.txt"
        rows = build(path, nan_rows=(2,))
        columns_tsv(cols)
        scales_txt(scales, empty=mode == "default")
        proc, pdf, bed = run_script(tmp_path, path, cols, scales, mode=mode)
        assert proc.returncode == 0, proc.stderr
        assert pdf.read_bytes().startswith(b"%PDF")
        got = bed_names(bed)
        order = expected_order(rows, GROUPS, range(4))
        assert [g[3] for g in got] == [rows[r][0][3] for r in order]
        assert [g[12] for g in got] == ["all_genes"] * 6 + ["all_genes_nostrand"] * 3

    def test_keep_and_theme(self, tmp_path):
        path, cols, scales = tmp_path / "m.gz", tmp_path / "cols.tsv", tmp_path / "s.txt"
        rows = build(path)
        columns_tsv(cols)
        scales_txt(scales)
        theme = tmp_path / "t.json"
        theme.write_text('{"theme": "presentation", "overrides": '
                         '{"heatmap_interpolation": {"mC": "bilinear"}, "profile_type": "se"}}')
        proc, pdf, bed = run_script(tmp_path, path, cols, scales, env="most", sort="no",
                                    theme=theme)
        assert proc.returncode == 0, proc.stderr
        assert [g[3] for g in bed_names(bed)] == [r[0][3] for r in rows]

    def test_bad_theme_stops(self, tmp_path):
        path, cols, scales = tmp_path / "m.gz", tmp_path / "cols.tsv", tmp_path / "s.txt"
        build(path)
        columns_tsv(cols)
        scales_txt(scales)
        theme = tmp_path / "t.json"
        theme.write_text('{"theme": "default", "overrides": {"heatmap_interpolation": {"mC": "cubic"}}}')
        proc, _, _ = run_script(tmp_path, path, cols, scales, theme=theme)
        assert proc.returncode != 0
        assert "heatmap_interpolation.mC" in proc.stderr

    def test_display_must_match_the_matrix(self, tmp_path):
        path, cols, scales = tmp_path / "m.gz", tmp_path / "cols.tsv", tmp_path / "s.txt"
        build(path)
        columns_tsv(cols, ALL_KEYS[:3])
        scales_txt(scales)
        proc, _, _ = run_script(tmp_path, path, cols, scales)
        assert proc.returncode != 0
        assert "Display rows do not match" in proc.stderr

    def test_images_are_drawn_at_pixel_size_without_smoothing(self, tmp_path):
        # One raster per column per region group, none interpolated by the viewer.
        path, cols, scales = tmp_path / "m.gz", tmp_path / "cols.tsv", tmp_path / "s.txt"
        build(path)
        columns_tsv(cols)
        scales_txt(scales)
        theme = tmp_path / "t.json"
        theme.write_text('{"theme": "default", "overrides": {"heatmap_height": 2, '
                         '"heatmap_column_width": 1, "heatmap_raster_dpi": 50}}')
        proc, pdf, _ = run_script(tmp_path, path, cols, scales, theme=theme)
        assert proc.returncode == 0, proc.stderr
        data = pdf.read_bytes()
        images = re.findall(rb"/Subtype /Image\s+/Width (\d+)\s+/Height (\d+)", data)
        heat = [(int(w), int(h)) for w, h in images if int(h) > 1]
        assert len(heat) == len(ALL_KEYS) * len(GROUPS)
        # 1 in x 50 dpi wide; the groups share the 2 in body less the gap
        # between them (0.08 in), in proportion to their regions (6 : 3).
        assert all(w == 50 for w, _ in heat)
        assert sorted({h for _, h in heat}) == [32, 64]
        assert b"/Interpolate true" not in data


@needs_ch
@needs_deeptools
def test_sorted_regions_match_plotheatmap(tmp_path):
    path, cols, scales = tmp_path / "m.gz", tmp_path / "cols.tsv", tmp_path / "s.txt"
    build(path, nan_rows=(1, 7))
    columns_tsv(cols)
    scales_txt(scales)
    proc, _, bed = run_script(tmp_path, path, cols, scales, env="most")
    assert proc.returncode == 0, proc.stderr
    dt_bed = tmp_path / "dt.bed"
    subprocess.run(["plotHeatmap", "-m", str(path), "-out", str(tmp_path / "dt.png"),
                    "--sortRegions", "descend", "--sortUsing", "mean",
                    "--outFileSortedRegions", str(dt_bed)], check=True, capture_output=True)
    ours, theirs = bed.read_text().splitlines(), dt_bed.read_text().splitlines()
    assert ours[0] == theirs[0]
    # Every column but blockStart, which deeptools writes as start minus the
    # first *digit* of the start (int() of the first character of a joined
    # string); it is 0 for a one-block region, as written here.
    cut = lambda lines: [ln.split("\t")[:11] + ln.split("\t")[12:] for ln in lines[1:]]  # noqa: E731
    assert cut(ours) == cut(theirs)
    assert all(ln.split("\t")[11] == "0" for ln in ours[1:])


# --- the rules ----------------------------------------------------------------

def _rule(name):
    text = SMK.read_text()
    block = text[text.index(f"rule {name}:"):]
    return block[:block.index("\nrule ", 1)]


def test_rule_draws_in_r_in_the_plots_env():
    block = _rule("plotting_heatmap_on_targetfile")
    assert "plotHeatmap" not in block and "plotheatmap_lowmem" not in block
    assert "conda: CONDA_ENV_PLOTS" in block
    assert "R_heatmap.R" in block
    assert "theme = PLOT_THEME_JSON" in block
    assert '"$TMPDIR/plot_theme.json"' in block
    assert 'display = lambda wildcards: define_key_for_plots(wildcards, "display")' in block
    assert "mem_mb = lambda wildcards, input, attempt: heatmap_mem_mb(input.matrix, attempt)" in block


def test_sorted_mc_heatmap_is_retired():
    text = SMK.read_text()
    for gone in ("rule sort_heatmap:", "rule plotting_sorted_heatmap_on_targetfile:",
                 "Heatmap_sorted__", "sorted_final_matrix_", "heatmap_sort_mc_after_others"):
        assert gone not in text


def test_default_targets_are_the_combined_heatmap():
    text = SMK.read_text()
    block = text[text.index("def define_final_combined_output"):]
    block = block[:block.index("\ndef ", 1)]
    assert "Heatmap__{mp}__all__" in block
    assert "Heatmap__regions__most__" not in block and "Heatmap__regions__mC__" not in block


def test_plots_env_has_complexheatmap():
    env = (REPO / "workflow" / "envs" / "epibutton_plots.yaml").read_text()
    assert "bioconductor-complexheatmap" in env
    assert "complexheatmap" not in (REPO / "workflow" / "envs" / "epibutton.yaml").read_text().lower()
