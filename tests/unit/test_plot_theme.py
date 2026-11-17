"""plot_theme.R: the shared theme every R plotting script takes its look from,
and plot_theme_spec, which hands the options-YAML choice to it as JSON.

The R tests source plot_theme.R in a throwaway Rscript and print what they
check. Skipped where Rscript with ggplot2 and jsonlite is not installed.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
THEME = REPO / "workflow" / "scripts" / "plot_theme.R"
MC_SCRIPT = REPO / "workflow" / "scripts" / "R_mc_profile_by_context.R"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"

sys.path.insert(0, os.path.join(REPO, "workflow", "scripts"))
from config_validation import plot_theme_spec  # noqa: E402


def _have_r():
    if not shutil.which("Rscript"):
        return False
    proc = subprocess.run(["Rscript", "-e", "library(ggplot2); library(jsonlite)"],
                          capture_output=True)
    return proc.returncode == 0


needs_r = pytest.mark.skipif(not _have_r(), reason="Rscript with ggplot2/jsonlite not found")


def _r(code, theme_spec=None, tmp_path=None):
    """Run R code after sourcing plot_theme.R; `th` is the loaded theme."""
    theme_file = ""
    if theme_spec is not None:
        theme_file = str(tmp_path / "theme.json")
        Path(theme_file).write_text(json.dumps(theme_spec))
    prog = f'source("{THEME}"); th <- epicc_theme("{theme_file}"); {code}'
    return subprocess.run(["Rscript", "-e", prog], capture_output=True, text=True)


# --- plot_theme_spec (Python side) -----------------------------------------

class TestPlotThemeSpec:
    def test_defaults_when_unset(self):
        assert json.loads(plot_theme_spec({})) == {"theme": "default", "overrides": {}}

    def test_passes_name_and_overrides(self):
        cfg = {"plot_theme": "presentation",
               "plot_theme_overrides": {"base_size": 12,
                                        "heatmap_colors": {"RNA": ["white", "darkgreen"]}}}
        assert json.loads(plot_theme_spec(cfg)) == {
            "theme": "presentation",
            "overrides": {"base_size": 12, "heatmap_colors": {"RNA": ["white", "darkgreen"]}}}

    def test_null_overrides_block_is_empty(self):
        # A YAML key with every child commented out parses as None.
        assert json.loads(plot_theme_spec({"plot_theme_overrides": None}))["overrides"] == {}

    def test_overrides_must_be_a_mapping(self):
        with pytest.raises(ValueError, match="plot_theme_overrides"):
            plot_theme_spec({"plot_theme_overrides": ["base_size", 12]})

    def test_name_must_be_a_string(self):
        with pytest.raises(ValueError, match="plot_theme"):
            plot_theme_spec({"plot_theme": 3})

    def test_stable_text_for_the_same_choice(self):
        # The JSON is a Snakemake param: key order must not make it differ
        # between runs and retrigger the plot rules.
        a = plot_theme_spec({"plot_theme_overrides": {"base_size": 12, "line_width": 1.0}})
        b = plot_theme_spec({"plot_theme_overrides": {"line_width": 1.0, "base_size": 12}})
        assert a == b


def test_profile_rule_hands_the_theme_to_the_mc_script():
    text = SMK.read_text()
    assert "theme = PLOT_THEME_JSON" in text
    assert '"$TMPDIR/plot_theme.json"' in text


# --- plot_theme.R ------------------------------------------------------------

@needs_r
class TestTheme:
    def test_default_theme_without_a_file(self, tmp_path):
        proc = _r("cat(th$base_size, th$line_width, length(th$palette))")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["11", "0.6", "20"]

    def test_palette_matches_the_snakemake_tab20(self, tmp_path):
        # Same list as _PALETTES["tab20"] in combined_analysis.smk.
        text = SMK.read_text()
        block = text[text.index('"tab20": ['):]
        block = block[:block.index("]")]
        smk_tab20 = [c.strip().strip('"').lower() for c in block.split("[", 1)[1].split(",") if c.strip()]
        proc = _r("cat(tolower(th$palette), sep='\\n')")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == smk_tab20

    def test_named_theme(self, tmp_path):
        proc = _r("cat(th$base_size)", {"theme": "presentation", "overrides": {}}, tmp_path)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.strip() == "15"

    def test_override_replaces_one_key(self, tmp_path):
        proc = _r("cat(th$base_size, th$line_width)",
                  {"theme": "default", "overrides": {"base_size": 14}}, tmp_path)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["14", "0.6"]

    def test_override_applies_on_top_of_the_named_theme(self, tmp_path):
        proc = _r("cat(th$base_size, th$line_width)",
                  {"theme": "presentation", "overrides": {"line_width": 2}}, tmp_path)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["15", "2"]

    def test_unknown_key_is_an_error(self, tmp_path):
        proc = _r("", {"theme": "default", "overrides": {"base_sise": 14}}, tmp_path)
        assert proc.returncode != 0
        assert "base_sise" in proc.stderr

    def test_wrong_type_is_an_error(self, tmp_path):
        proc = _r("", {"theme": "default", "overrides": {"base_size": "big"}}, tmp_path)
        assert proc.returncode != 0
        assert "base_size" in proc.stderr

    def test_unknown_theme_is_an_error(self, tmp_path):
        proc = _r("", {"theme": "neon", "overrides": {}}, tmp_path)
        assert proc.returncode != 0
        assert "neon" in proc.stderr and "presentation" in proc.stderr

    def test_heatmap_colors_is_an_open_map(self, tmp_path):
        # New data types may be added; others keep their defaults.
        spec = {"theme": "default",
                "overrides": {"heatmap_colors": {"H3K9me2": ["white", "black"]}}}
        proc = _r("cat(epicc_heatmap_colors(th, 'H3K9me2'), '|', "
                  "epicc_heatmap_colors(th, 'mC')[1])", spec, tmp_path)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["white", "black", "|", "#FFF5EB"]

    def test_heatmap_colour_needs_two_colours(self, tmp_path):
        spec = {"theme": "default", "overrides": {"heatmap_colors": {"RNA": "red"}}}
        proc = _r("", spec, tmp_path)
        assert proc.returncode != 0
        assert "heatmap_colors.RNA" in proc.stderr

    @pytest.mark.parametrize("dtype,first", [
        ("mCHG", "#FFF5EB"),        # context -> mC
        ("ChIP_broad", "#00004C"),  # combined token -> its family
        ("H3K27me3", "#00004C"),    # unknown -> default
    ])
    def test_heatmap_colour_fallbacks(self, tmp_path, dtype, first):
        proc = _r(f"cat(epicc_heatmap_colors(th, '{dtype}')[1])")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.strip() == first

    def test_dashes_split_only_shared_colours(self, tmp_path):
        proc = _r("d <- epicc_dashes(th, c(WT='#1f77b4', ddm1='#ff7f0e', WT_CR='#1f77b4')); "
                  "cat(paste(names(d), d, sep='='))")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["WT=solid", "ddm1=solid", "WT_CR=dashed"]

    def test_legend_wraps_at_max_rows(self, tmp_path):
        proc = _r("cat(epicc_legend_ncol(th, 1), epicc_legend_ncol(th, 25), epicc_legend_ncol(th, 26))")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["1", "1", "2"]

    def test_ggtheme_builds(self, tmp_path):
        proc = _r("g <- epicc_ggtheme(th); cat(inherits(g, 'theme'), "
                  "inherits(g$panel.grid.minor, 'element_blank'))")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["TRUE", "TRUE"]


@needs_r
def test_mc_profile_takes_a_theme_file(tmp_path):
    ticks = ["-2.0Kb", "", "start", "end", "", "2.0Kb"]
    labels = ["WT_mCG", "ddm1_mCG"]
    lines = ["bin labels\t\t" + "\t".join(ticks * len(labels)),
             "bins\t\t" + "\t".join(str(float(i + 1)) for i in range(6))]
    lines += [f"{l}\tall_genes\t" + "\t".join(str(10.0 + i) for i in range(6)) for l in labels]
    data, samples, theme, out = (tmp_path / n for n in ("p.tab", "s.tsv", "t.json", "out.pdf"))
    data.write_text("\n".join(lines) + "\n")
    samples.write_text("".join(f"{l}\tmCG\t#1f77b4\t{l.split('_')[0]}\n" for l in labels))

    theme.write_text(json.dumps({"theme": "presentation", "overrides": {"line_width": 1.5}}))
    proc = subprocess.run(["Rscript", str(MC_SCRIPT), str(data), str(samples), str(out), "y", str(theme)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert out.read_bytes().startswith(b"%PDF")

    theme.write_text(json.dumps({"theme": "default", "overrides": {"nope": 1}}))
    proc = subprocess.run(["Rscript", str(MC_SCRIPT), str(data), str(samples), str(out), "y", str(theme)],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert "nope" in proc.stderr
