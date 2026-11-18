"""R_profile.R / plot_profile.R: Profile__ and Profile_pergroup__ drawn in R.

The matrices are written in deeptools' format (test_matrix_io.write_matrix):
a `most`-style matrix with interleaved stranded columns and a nostrand group,
and an mC matrix (sample x context columns, NaNs). The checks look inside the
built ggplots (ggplot_build) rather than at pixels: which y range each panel
got, and which values each line draws.
"""

import json
import math
import re
import subprocess
from pathlib import Path

import pytest

from .test_matrix_io import needs_r, write_matrix

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"
SCRIPT = SCRIPTS / "R_profile.R"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"
NAN = float("nan")
NBINS = 4

MOST_KEYS = ["WT__H3K9me2", "ddm1__H3K9me2", "WT__RNA_plus", "WT__RNA_minus",
             "ddm1__RNA_plus", "ddm1__RNA_minus", "WT__H3K4me3"]
MOST_GROUPS = [("genes", 5), ("genes_nostrand", 3)]
MC_KEYS = ["WT__mCG", "ddm1__mCG", "WT__mCHG", "ddm1__mCHG"]
COLOURS = {"WT": "#1f77b4", "ddm1": "#ff7f0e"}


def _level(key):
    mark = key.split("__")[-1]
    base = {"H3K9me2": 50, "RNA_plus": 8, "RNA_minus": 1, "H3K4me3": 0.05,
            "mCG": 60, "mCHG": 20}[mark]
    return base * (0.5 if key.startswith("ddm1") else 1)


def build(path, keys, groups, nan_every=0, nbins=NBINS, bin_size=100):
    rows, r = [], 0
    for g, n in groups:
        for i in range(n):
            vals = []
            for k in keys:
                lv = _level(k)
                v = [lv * (1 + ((r * 3 + b * 5 + len(k)) % 7) / 7) for b in range(nbins)]
                if nan_every and r % nan_every == 0:
                    v[1] = NAN
                vals += v
            strand = "." if g.endswith("nostrand") else "+-"[i % 2]
            rows.append((("chr1", r * 100, r * 100 + 50, f"{g}{i}", 0, strand), vals))
            r += 1
    write_matrix(path, groups, keys, rows, body=nbins * bin_size, bin_size=bin_size,
                 upstream=0, downstream=0)
    return rows


def columns_tsv(path, keys):
    lines = []
    for k in keys:
        sample, track = k.split("__")
        lines.append(f"{k}\t{sample} {track}\t{COLOURS[sample]}\t{track}\t{sample}")
    path.write_text("\n".join(lines) + "\n")


def run_r(code):
    prog = "".join(f'source("{SCRIPTS / f}")\n' for f in
                   ("matrix_io.R", "plot_scales.R", "plot_theme.R", "plot_profile.R")) + code
    return subprocess.run(["Rscript", "-e", prog], capture_output=True, text=True)


def run_r_ok(code):
    proc = run_r(code)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


def run_script(tmp, matrix, cols, env="most", mode="type", average="mean", theme=None):
    out1, out2 = tmp / "p1.pdf", tmp / "p2.pdf"
    args = ["Rscript", str(SCRIPT), str(matrix), str(cols), "regions", env, mode, average,
            str(out1), str(out2)]
    if theme is not None:
        args.append(str(theme))
    return subprocess.run(args, capture_output=True, text=True), out1, out2


@pytest.fixture
def most(tmp_path):
    path, cols = tmp_path / "most.gz", tmp_path / "cols.tsv"
    rows = build(path, MOST_KEYS, MOST_GROUPS)
    columns_tsv(cols, MOST_KEYS)
    return path, cols, rows


@pytest.fixture
def mc(tmp_path):
    path, cols = tmp_path / "mc.gz", tmp_path / "mc_cols.tsv"
    rows = build(path, MC_KEYS, [("all_TEs", 4), ("all_TEs_nostrand", 2)], nan_every=3)
    columns_tsv(cols, MC_KEYS)
    return path, cols, rows


# Every panel's y range and its facet value, from ggplot_build.
PANEL_RANGES = '''
panel_ranges <- function(p, var) {
  b <- ggplot_build(p)
  lay <- b$layout$layout
  for (k in seq_len(nrow(lay))) {
    r <- b$layout$panel_params[[k]]$y.range
    cat(as.character(lay[[var]][k]), sprintf("%.10g", r), "\\n", sep = "\\t")
  }
}
display <- function(path) read.delim(path, header = FALSE, sep = "\\t", quote = "",
  col.names = c("key", "label", "colour", "track", "series"), colClasses = "character")
'''


def parse_ranges(out):
    res = {}
    for line in out.strip().splitlines():
        f = [x for x in line.split("\t") if x.strip()]
        res.setdefault(f[0], []).append((float(f[1]), float(f[2])))
    return res


@needs_r
class TestScript:
    @pytest.mark.parametrize("mode", ["type", "sample", "default"])
    @pytest.mark.parametrize("average", ["mean", "median"])
    def test_draws_both_pdfs(self, most, tmp_path, mode, average):
        path, cols, _ = most
        proc, p1, p2 = run_script(tmp_path, path, cols, mode=mode, average=average)
        assert proc.returncode == 0, proc.stderr
        assert p1.read_bytes().startswith(b"%PDF") and p2.read_bytes().startswith(b"%PDF")

    @pytest.mark.parametrize("mode", ["type", "default"])
    def test_mc(self, mc, tmp_path, mode):
        path, cols, _ = mc
        proc, p1, p2 = run_script(tmp_path, path, cols, env="mC", mode=mode)
        assert proc.returncode == 0, proc.stderr
        assert p2.read_bytes().startswith(b"%PDF")

    def test_one_context_whose_labels_dropped_it(self, tmp_path):
        # methylation_contexts: [CG]: the series is the whole label.
        path, cols = tmp_path / "cg.gz", tmp_path / "cg.tsv"
        build(path, ["WT__mCG", "ddm1__mCG"], [("all_genes", 3)])
        cols.write_text("WT__mCG\tWT\t#1f77b4\tmCG\tWT\nddm1__mCG\tddm1\t#ff7f0e\tmCG\tddm1\n")
        proc, _, p2 = run_script(tmp_path, path, cols, env="mC")
        assert proc.returncode == 0, proc.stderr

    def test_display_must_match_the_matrix(self, most, tmp_path):
        path, cols, _ = most
        columns_tsv(cols, MOST_KEYS[:-1])
        proc, _, _ = run_script(tmp_path, path, cols)
        assert proc.returncode != 0
        assert "do not match the matrix columns" in proc.stderr

    def test_mc_needs_a_series(self, mc, tmp_path):
        path, cols, _ = mc
        cols.write_text("".join(f"{k}\t{k}\t#1f77b4\t{k.split('__')[1]}\t\n" for k in MC_KEYS))
        proc, _, _ = run_script(tmp_path, path, cols, env="mC")
        assert proc.returncode != 0
        assert "series" in proc.stderr

    def test_theme_file(self, most, tmp_path):
        path, cols, _ = most
        theme = tmp_path / "t.json"
        theme.write_text(json.dumps({"theme": "presentation",
                                     "overrides": {"line_width": 1.5, "profile_type": "se"}}))
        proc, _, _ = run_script(tmp_path, path, cols, theme=theme)
        assert proc.returncode == 0, proc.stderr
        theme.write_text(json.dumps({"theme": "default", "overrides": {"nope": 1}}))
        proc, _, _ = run_script(tmp_path, path, cols, theme=theme)
        assert proc.returncode != 0
        assert "nope" in proc.stderr

    @pytest.mark.parametrize("ptype", ["fill", "se", "std"])
    def test_profile_types(self, most, tmp_path, ptype):
        path, cols, _ = most
        theme = tmp_path / "t.json"
        theme.write_text(json.dumps({"theme": "default", "overrides": {"profile_type": ptype}}))
        proc, p1, _ = run_script(tmp_path, path, cols, theme=theme)
        assert proc.returncode == 0, proc.stderr


@needs_r
class TestLimits:
    @pytest.mark.parametrize("mode", ["type", "sample"])
    def test_each_column_panel_gets_its_own_columns_limits(self, most, mode):
        """Profile__ on an interleaved stranded matrix: panel `key` spans
        exactly scale_for_column(key)'s y range -- looked up by key, not by
        position within a mark."""
        path, cols, _ = most
        out = run_r_ok(PANEL_RANGES + f'''
            m <- read_matrix("{path}")
            lim <- profile_limits(m, profile_data(m, "mean"), "{mode}")
            s <- matrix_scales(m, "{mode}", "mean")
            p <- plot_profile_per_column(m, display("{cols}"), epicc_theme(), lim)$plot
            panel_ranges(p, "panel")
            for (k in m$samples$label) {{
              r <- scale_for_column(s, k)
              cat("want", k, sprintf("%.10g", c(r$ymin, r$ymax)), "\\n", sep = "\\t")
            }}''')
        got, want = {}, {}
        for line in out.strip().splitlines():
            f = [x for x in line.split("\t") if x.strip()]
            if f[0] == "want":
                want[f[1]] = (float(f[2]), float(f[3]))
            else:
                got[f[0]] = (float(f[1]), float(f[2]))
        assert list(got) == MOST_KEYS
        for k in MOST_KEYS:
            assert got[k] == pytest.approx(want[k], rel=1e-9), k
        # type mode: RNA_plus and RNA_minus are different marks
        if mode == "type":
            assert got["WT__RNA_plus"] == got["ddm1__RNA_plus"] != got["WT__RNA_minus"]

    def test_per_group_shares_one_range(self, most):
        path, cols, _ = most
        out = run_r_ok(PANEL_RANGES + f'''
            m <- read_matrix("{path}")
            lim <- profile_limits(m, profile_data(m, "mean"), "type")
            p <- plot_profile_per_group(m, display("{cols}"), epicc_theme(), lim)$plot
            panel_ranges(p, "region")
            cat("want", sprintf("%.10g", shared_ylim(lim)), "\\n", sep = "\\t")''')
        r = parse_ranges(out)
        want = r.pop("want")[0]
        assert list(r) == ["genes(5)", "genes_nostrand(3)"]
        for rng in r.values():
            assert rng[0] == pytest.approx(want, rel=1e-9)

    def test_default_mode_is_one_shared_range_not_99999(self, most):
        """heatmaps_scales: default once collapsed to --yMin 99999 --yMax
        -99999; now every panel shares default_scales()' range."""
        path, cols, _ = most
        out = run_r_ok(PANEL_RANGES + f'''
            m <- read_matrix("{path}")
            lim <- profile_limits(m, profile_data(m, "mean"), "default")
            d <- default_scales(m, "mean")
            p <- plot_profile_per_column(m, display("{cols}"), epicc_theme(), lim)$plot
            panel_ranges(p, "panel")
            cat("want", sprintf("%.10g", c(d$ymin, d$ymax)), "\\n", sep = "\\t")''')
        r = parse_ranges(out)
        want = r.pop("want")[0]
        assert len({v[0] for v in r.values()}) == 1
        assert next(iter(r.values()))[0] == pytest.approx(want, rel=1e-9)
        assert abs(want[0]) < 1000 and abs(want[1]) < 1000

    def test_mc_context_panels_share_their_contexts_range(self, mc):
        path, cols, _ = mc
        out = run_r_ok(PANEL_RANGES + f'''
            m <- read_matrix("{path}")
            lim <- profile_limits(m, profile_data(m, "mean"), "sample")
            dsp <- display("{cols}")
            p <- plot_profile_by_context(m, dsp, epicc_theme(), lim)$plot
            panel_ranges(p, "context")
            for (ctx in c("mCG", "mCHG"))
              cat("want", ctx, sprintf("%.10g", shared_ylim(lim[dsp$track == ctx, ])), "\\n", sep = "\\t")''')
        got, want = {}, {}
        for line in out.strip().splitlines():
            f = [x for x in line.split("\t") if x.strip()]
            if f[0] == "want":
                want[f[1]] = (float(f[2]), float(f[3]))
            else:
                got.setdefault(f[0], []).append((float(f[1]), float(f[2])))
        assert set(got) == {"mCG", "mCHG"}
        for ctx, ranges in got.items():
            assert len(ranges) == 2  # one panel per region group
            for rng in ranges:
                assert rng == pytest.approx(want[ctx], rel=1e-9)

    def test_default_mode_profile_panels_each_autoscale(self, mc):
        """Profile__ under default: each panel on its own curves' range, so
        mCHG isn't drawn on mCG's (the per-context panels autoscale too)."""
        path, cols, _ = mc
        out = run_r_ok(PANEL_RANGES + f'''
            m <- read_matrix("{path}")
            d <- profile_data(m, "mean")
            lim <- profile_limits(m, d, "default", per_column = TRUE)
            p <- plot_profile_per_column(m, display("{cols}"), epicc_theme(), lim, data = d)$plot
            panel_ranges(p, "panel")
            for (k in m$samples$label)
              cat("want", k, sprintf("%.10g", default_ylim(d[d$sample == k, ])), "\\n", sep = "\\t")''')
        got, want = {}, {}
        for line in out.strip().splitlines():
            f = [x for x in line.split("\t") if x.strip()]
            if f[0] == "want":
                want[f[1]] = (float(f[2]), float(f[3]))
            else:
                got[f[0]] = (float(f[1]), float(f[2]))
        for k in MC_KEYS:
            assert got[k] == pytest.approx(want[k], rel=1e-9), k
        assert got["WT__mCHG"][1] < got["WT__mCG"][0]


@needs_r
class TestAllEnv:
    """The `all` matrix is matrix_join_all.py's join of a `most` and an mC
    matrix: repeated bins at the common bin size, 0-filled non-mC rows and
    NaN mC values side by side."""

    KEYS = ["WT__H3K9me2", "WT__H3K4me3"] + MC_KEYS

    @pytest.fixture
    def joined(self, tmp_path):
        groups = [("genes", 5), ("genes_nostrand", 3)]
        most, mc = tmp_path / "most.gz", tmp_path / "mc.gz"
        # most at bin size 50, mC at 100: the join repeats each mC bin twice.
        build(most, self.KEYS[:2], groups, nbins=2 * NBINS, bin_size=50)
        build(mc, MC_KEYS, groups, nan_every=3)
        out, cols = tmp_path / "all.gz", tmp_path / "all_cols.tsv"
        subprocess.run(["python3", str(SCRIPTS / "matrix_join_all.py"), "-o", str(out),
                        "--most", str(most), "--mc", str(mc)], check=True)
        columns_tsv(cols, self.KEYS)
        return out, cols

    @pytest.mark.parametrize("mode", ["type", "sample", "default"])
    def test_script_draws_the_joined_matrix(self, tmp_path, joined, mode):
        path, cols = joined
        proc, out1, out2 = run_script(tmp_path, path, cols, env="all", mode=mode)
        assert proc.returncode == 0, proc.stderr
        assert out1.stat().st_size > 0 and out2.stat().st_size > 0

    def test_pergroup_has_a_panel_per_mark(self, joined):
        """mC and ChIP each keep their own range in Profile_pergroup__."""
        path, cols = joined
        out = run_r_ok(PANEL_RANGES + f'''
            m <- read_matrix("{path}")
            d <- profile_data(m, "mean")
            p <- plot_profile_by_context(m, display("{cols}"), epicc_theme(),
                                         profile_limits(m, d, "type"), data = d)$plot
            panel_ranges(p, "context")''')
        r = parse_ranges(out)
        assert set(r) == {"H3K9me2", "H3K4me3", "mCG", "mCHG"}
        assert all(len(v) == 2 for v in r.values())  # one per region group
        assert all(hi < 1 for _, hi in r["H3K4me3"])
        assert all(lo > 20 for lo, _ in r["mCG"])


@needs_r
class TestClipping:
    """deeptools applied --yMin/--yMax after drawing and clipped. A band (se,
    std) or a fill baseline must not widen a panel past its limits: two
    panels of one mark keep one axis, and an mC fill is not dragged to 0."""

    @pytest.mark.parametrize("ptype", ["se", "std", "fill"])
    def test_type_mode_panels_keep_their_marks_range(self, most, ptype):
        path, cols, _ = most
        out = run_r_ok(PANEL_RANGES + f'''
            th <- epicc_theme(); th$profile_type <- "{ptype}"
            m <- read_matrix("{path}")
            d <- profile_data(m, "mean", .spread(th))
            lim <- profile_limits(m, d, "type")
            s <- matrix_scales(m, "type", "mean")
            p <- plot_profile_per_column(m, display("{cols}"), th, lim, data = d)$plot
            panel_ranges(p, "panel")
            for (k in m$samples$label) {{
              r <- scale_for_column(s, k)
              cat("want", k, sprintf("%.10g", c(r$ymin, r$ymax)), "\\n", sep = "\\t")
            }}''')
        got, want = {}, {}
        for line in out.strip().splitlines():
            f = [x for x in line.split("\t") if x.strip()]
            if f[0] == "want":
                want[f[1]] = (float(f[2]), float(f[3]))
            else:
                got[f[0]] = (float(f[1]), float(f[2]))
        for k in MOST_KEYS:
            assert got[k] == pytest.approx(want[k], rel=1e-9), k
        assert got["WT__H3K9me2"] == got["ddm1__H3K9me2"]

    def test_band_is_clipped_to_the_panel(self, most):
        path, cols, _ = most
        out = run_r_ok(f'''
            th <- epicc_theme(); th$profile_type <- "std"
            m <- read_matrix("{path}")
            d <- profile_data(m, "mean", "std")
            lim <- profile_limits(m, d, "type")
            p <- plot_profile_per_column(m, read.delim("{cols}", header = FALSE,
                   col.names = c("key", "label", "colour", "track", "series")), th, lim, data = d)$plot
            r <- layer_data(p, 1)
            b <- ggplot_build(p)$layout$layout
            k <- as.character(b$panel[match(r$PANEL, b$PANEL)])
            l <- lim[match(k, lim$key), ]
            cat(all(r$ymin >= l$ymin - 1e-12 & r$ymax <= l$ymax + 1e-12),
                any(d$lo < lim$ymin[match(d$sample, lim$key)]))''')
        # the unclipped std band does reach below the limits; drawn, it doesn't
        assert out.split() == ["TRUE", "TRUE"]

    @pytest.mark.parametrize("mode", ["type", "default"])
    def test_mc_fill_is_not_dragged_to_zero(self, mc, mode):
        path, cols, _ = mc
        out = run_r_ok(PANEL_RANGES + f'''
            th <- epicc_theme(); th$profile_type <- "fill"
            m <- read_matrix("{path}")
            d <- profile_data(m, "mean")
            lim <- profile_limits(m, d, "{mode}")
            p <- plot_profile_by_context(m, display("{cols}"), th, lim, use_limits = "{mode}" != "default",
                                         data = d)$plot
            panel_ranges(p, "context")''')
        r = parse_ranges(out)
        # mCG sits at 30-120: its panels must not start at 0
        assert all(lo > 20 for lo, _ in r["mCG"])


@needs_r
class TestCurves:
    @pytest.mark.parametrize("average", ["mean", "median"])
    def test_lines_are_profiles_scale(self, most, average):
        """Profile__ draws profiles_scale (deeptools drew it as a mean even
        under median, while its limits came from the median)."""
        path, cols, rows = most
        out = run_r_ok(PANEL_RANGES + f'''
            m <- read_matrix("{path}")
            lim <- profile_limits(m, profile_data(m, "{average}"), "type")
            p <- plot_profile_per_column(m, display("{cols}"), epicc_theme(), lim, "{average}")$plot
            d <- layer_data(p, 1)
            b <- ggplot_build(p)$layout$layout
            k <- b$panel[match(d$PANEL, b$PANEL)]
            sel <- k == "ddm1__RNA_minus" & d$group == min(d$group[k == "ddm1__RNA_minus"])
            cat(sprintf("%.10g", d$y[sel][order(d$x[sel])]))''')
        col = MOST_KEYS.index("ddm1__RNA_minus")
        n = MOST_GROUPS[0][1]
        want = []
        for b in range(NBINS):
            v = sorted(rows[r][1][col * NBINS + b] for r in range(n))
            want.append(sum(v) / n if average == "mean" else
                        (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2))
        assert [float(x) for x in out.split()] == pytest.approx(want, rel=1e-6)

    def test_se_band(self, mc):
        """se = nanstd (population) / sqrt(regions with a value), per bin."""
        path, _, rows = mc
        out = run_r_ok(f'''
            m <- read_matrix("{path}")
            d <- profile_data(m, "mean", "se")
            s <- d[d$sample == "WT__mCG" & d$group == "all_TEs", ]
            cat(sprintf("%.10g", s$hi - s$value))''')
        got = [float(x) for x in out.split()]
        for b in range(NBINS):
            v = [rows[r][1][b] for r in range(4)]
            v = [x for x in v if not math.isnan(x)]
            mu = sum(v) / len(v)
            sd = math.sqrt(sum((x - mu) ** 2 for x in v) / len(v))
            assert got[b] == pytest.approx(sd / math.sqrt(len(v)), rel=1e-6)


def test_rule_draws_in_r_in_the_plots_env():
    text = SMK.read_text()
    block = text[text.index("rule plotting_profile_on_targetfile:"):]
    block = block[:block.index("\nrule ", 1)]
    assert "plotProfile" not in block
    assert "conda: CONDA_ENV_PLOTS" in block
    assert "R_profile.R" in block
    # Labels and colours stay out of every data rule; the plot rule takes them.
    assert re.search(r"display = lambda wildcards: define_key_for_plots\(wildcards, \"display\"\)",
                     block)
