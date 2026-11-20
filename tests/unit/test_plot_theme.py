"""plot_theme.R: the shared theme every R plotting script takes its look from,
and plot_theme_spec, which hands the options-YAML choice to it as JSON.

The R tests source plot_theme.R in a throwaway Rscript and print what they
check. Skipped where Rscript with ggplot2 and jsonlite is not installed.
"""

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
THEME = REPO / "workflow" / "scripts" / "plot_theme.R"
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"

sys.path.insert(0, os.path.join(REPO, "workflow", "scripts"))
from config_validation import legacy_plot_params, plot_theme_spec  # noqa: E402


def spec_json(cfg, r_drawn=()):
    text, _ = plot_theme_spec(cfg, r_drawn)
    return json.loads(text)


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
        assert spec_json({}) == {"theme": "default", "overrides": {}}

    def test_passes_name_and_overrides(self):
        cfg = {"plot_theme": "presentation",
               "plot_theme_overrides": {"base_size": 12,
                                        "heatmap_colors": {"RNA": ["white", "darkgreen"]}}}
        assert spec_json(cfg) == {
            "theme": "presentation",
            "overrides": {"base_size": 12, "heatmap_colors": {"RNA": ["white", "darkgreen"]}}}

    def test_null_overrides_block_is_empty(self):
        # A YAML key with every child commented out parses as None.
        assert spec_json({"plot_theme_overrides": None})["overrides"] == {}

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


# --- deprecated heatmaps_plot_params / profiles_plot_params -------------------
# (the heatmaps half restored from branch claude/r-plot-legacy-flags, fc1130c,
# with --interpolationMethod now translated rather than dropped)

# What the options file shipped before both families moved to R.
SHIPPED_HEATMAP_PARAMS = {
    env: "--colorMap seismic --interpolationMethod bilinear"
    for env in ("all", "most", "ChIP", "TF", "RNA", "sRNA", "ATAC")
}
SHIPPED_HEATMAP_PARAMS["mC"] = "--colorMap Oranges --missingDataColor 0.9 --interpolationMethod nearest"
SHIPPED = {"heatmaps_plot_params": SHIPPED_HEATMAP_PARAMS,
           "profiles_plot_params": "--plotType lines"}
PROFILES = ("profiles",)
HEATMAPS = ("heatmaps",)
BOTH = ("heatmaps", "profiles")


class TestLegacyPlotParams:
    def test_shipped_options_file(self):
        import yaml
        cfg = yaml.safe_load((REPO / "config" / "epicc-options.yaml").read_text())
        # Both families are drawn in R: the keys are gone from the shipped
        # file and only translated when an older options file still has them.
        for key in ("profiles_plot_params", "heatmaps_plot_params",
                    "heatmap_sort_mc_after_others"):
            assert key not in cfg

    def test_nothing_translated_while_deeptools_draws(self):
        cfg = {"heatmaps_plot_params": {"RNA": "--colorMap viridis --kmeans 3"},
               "profiles_plot_params": "--plotType overlapped_lines",
               "heatmap_sort_mc_after_others": True}
        assert legacy_plot_params(cfg) == ({}, [])
        assert spec_json(cfg)["overrides"] == {}

    def test_old_defaults_add_no_overrides(self):
        overrides, warnings = legacy_plot_params(SHIPPED, BOTH)
        assert overrides == {}
        assert len(warnings) == 2
        assert all("deprecated" in w and "match the theme defaults" in w for w in warnings)

    def test_colour_maps_lists_and_missing_colour(self):
        cfg = {"heatmaps_plot_params": {
            "RNA": "--colorMap viridis --interpolationMethod bilinear",
            "ChIP": "--colorList white,darkgreen",
            "TF": "--colorMap seismic",
            "mC": "--colorMap Reds_r --missingDataColor 0.5"}}
        overrides, warnings = legacy_plot_params(cfg, HEATMAPS)
        assert overrides == {
            "heatmap_colors": {"RNA": "viridis", "ChIP": ["white", "darkgreen"], "mC": "Reds_r"},
            "heatmap_missing_color": "#808080"}
        assert "Translated to plot_theme_overrides" in warnings[0]

    def test_all_and_most_set_every_non_mc_type(self):
        # The theme has its own entry per non-mC type, which the lookup tries
        # before 'default': 'all'/'most' setting the fallback alone would
        # change no column.
        cfg = {"heatmaps_plot_params": {"all": "--colorMap viridis --interpolationMethod nearest",
                                        "most": "--colorMap=viridis"}}
        overrides, _ = legacy_plot_params(cfg, HEATMAPS)
        non_mc = ("default", "ChIP", "TF", "ATAC", "RNA", "sRNA")
        assert overrides == {"heatmap_colors": {t: "viridis" for t in non_mc},
                             "heatmap_interpolation": {t: "nearest" for t in non_mc}}
        cfg = {"heatmaps_plot_params": {"all": "--colorMap viridis", "most": "--colorMap magma"}}
        with pytest.raises(ValueError, match="both set heatmap_colors of every non-mC type"):
            legacy_plot_params(cfg, HEATMAPS)

    def test_an_envs_own_entry_wins_over_most(self):
        cfg = {"heatmaps_plot_params": {"most": "--colorMap viridis", "RNA": "--colorMap magma",
                                        "mC": "--colorMap Oranges"}}
        overrides, _ = legacy_plot_params(cfg, HEATMAPS)
        assert overrides["heatmap_colors"]["RNA"] == "magma"
        assert overrides["heatmap_colors"]["ChIP"] == "viridis"
        assert "mC" not in overrides["heatmap_colors"]

    def test_non_mc_types_match_the_theme(self):
        # config_validation's list must name every non-mC entry plot_theme.R
        # gives heatmap_colors and heatmap_interpolation.
        from config_validation import _NON_MC_TYPES
        text = THEME.read_text()
        for key in ("heatmap_colors", "heatmap_interpolation"):
            block = text[text.index(f"  {key} = list(") + len(f"  {key} = list("):]
            block = block[:block.index("\n  )")]
            entries = re.findall(r"^\s+(\w+) =", block, flags=re.M)
            assert sorted(e for e in entries if e != "mC") == sorted(_NON_MC_TYPES), key

    @needs_r
    def test_most_flags_reach_the_chip_columns(self, tmp_path):
        text, _ = plot_theme_spec({"heatmaps_plot_params": {"most": "--interpolationMethod nearest"}},
                                  HEATMAPS)
        out = _r('cat(epicc_heatmap_interpolation(th, "ChIP"), '
                 'epicc_heatmap_interpolation(th, "mCG"))', json.loads(text), tmp_path)
        assert out.stdout.split() == ["nearest", "nearest"]

    def test_interpolation(self):
        cfg = {"heatmaps_plot_params": {"ChIP": "--interpolationMethod nearest",
                                        "mC": "--interpolationMethod bilinear",
                                        "RNA": "--interpolationMethod bilinear"}}
        overrides, _ = legacy_plot_params(cfg, HEATMAPS)
        # RNA bilinear is the default and is dropped.
        assert overrides == {"heatmap_interpolation": {"ChIP": "nearest", "mC": "bilinear"}}

    @pytest.mark.parametrize("method,mode", [("none", "nearest"), ("bicubic", "bilinear"),
                                             ("lanczos", "bilinear"), ("antialiased", "bilinear")])
    def test_other_matplotlib_methods_map_with_a_note(self, method, mode):
        cfg = {"heatmaps_plot_params": {"sRNA": f"--interpolationMethod {method}"}}
        overrides, warnings = legacy_plot_params(cfg, HEATMAPS)
        assert overrides.get("heatmap_interpolation", {}).get("sRNA", "bilinear") == mode
        assert f"--interpolationMethod {method} is drawn as '{mode}'" in warnings[0]

    def test_unknown_interpolation_method(self):
        with pytest.raises(ValueError, match="not a matplotlib interpolation method"):
            legacy_plot_params({"heatmaps_plot_params": {"RNA": "--interpolationMethod smooth"}},
                               HEATMAPS)

    def test_unknown_colour_map(self):
        with pytest.raises(ValueError, match="unknown colour map 'nonesuch'"):
            legacy_plot_params({"heatmaps_plot_params": {"RNA": "--colorMap nonesuch"}}, HEATMAPS)

    def test_flag_without_theme_equivalent_stops_the_run(self):
        with pytest.raises(ValueError, match=r"heatmaps_plot_params\.RNA: '--kmeans'"):
            legacy_plot_params({"heatmaps_plot_params": {"RNA": "--colorMap seismic --kmeans 3"}},
                               HEATMAPS)
        with pytest.raises(ValueError, match="'--yAxisLabel'"):
            legacy_plot_params({"profiles_plot_params": "--yAxisLabel x"}, PROFILES)

    def test_per_sample_colour_maps_are_refused(self):
        with pytest.raises(ValueError, match="one colour map per data type"):
            legacy_plot_params({"heatmaps_plot_params": {"RNA": "--colorMap Reds Blues"}},
                               HEATMAPS)

    def test_conflicting_missing_colours(self):
        cfg = {"heatmaps_plot_params": {"RNA": "--missingDataColor black",
                                        "mC": "--missingDataColor white"}}
        with pytest.raises(ValueError, match="differs between envs"):
            legacy_plot_params(cfg, HEATMAPS)

    def test_heatmap_params_must_be_a_mapping(self):
        with pytest.raises(ValueError, match="must map each env"):
            legacy_plot_params({"heatmaps_plot_params": "--colorMap viridis"}, HEATMAPS)

    def test_sort_mc_after_others_is_reported_as_ignored(self):
        for value in (True, False):
            overrides, warnings = legacy_plot_params(
                {"heatmap_sort_mc_after_others": value}, HEATMAPS)
            assert overrides == {}
            assert len(warnings) == 1
            assert "heatmap_sort_mc_after_others is deprecated and ignored" in warnings[0]
            assert "--env mC" in warnings[0]

    def test_plot_type(self):
        overrides, warnings = legacy_plot_params({"profiles_plot_params": "--plotType se"}, PROFILES)
        assert overrides == {"profile_type": "se"}
        assert 'profile_type="se"' in warnings[0]
        with pytest.raises(ValueError, match="--plotType heatmap is not drawn in R"):
            legacy_plot_params({"profiles_plot_params": "--plotType heatmap"}, PROFILES)

    @pytest.mark.parametrize("flags", ["--plotType=se", "--plotType se", "--plotType=se "])
    def test_plot_type_equals_form(self, flags):
        # argparse (and so deeptools) took --flag=value as well as --flag value
        assert legacy_plot_params({"profiles_plot_params": flags}, PROFILES)[0] == \
            {"profile_type": "se"}

    def test_equals_form_of_an_unknown_flag_names_the_flag(self):
        with pytest.raises(ValueError, match="'--yAxisLabel'"):
            legacy_plot_params({"profiles_plot_params": "--yAxisLabel=x"}, PROFILES)

    def test_only_drawn_families_are_translated(self):
        cfg = {"heatmaps_plot_params": {"RNA": "--colorMap viridis"},
               "profiles_plot_params": "--plotType fill"}
        assert legacy_plot_params(cfg, PROFILES)[0] == {"profile_type": "fill"}
        assert legacy_plot_params(cfg, HEATMAPS)[0] == {"heatmap_colors": {"RNA": "viridis"}}

    def test_user_overrides_win(self):
        cfg = {"heatmaps_plot_params": {"RNA": "--colorMap viridis", "sRNA": "--colorMap magma",
                                        "ChIP": "--interpolationMethod nearest"},
               "profiles_plot_params": "--plotType fill",
               "plot_theme_overrides": {"profile_type": "se",
                                        "heatmap_colors": {"RNA": ["white", "red"]},
                                        "heatmap_interpolation": {"mC": "bilinear"}}}
        assert spec_json(cfg, BOTH)["overrides"] == {
            "profile_type": "se",
            "heatmap_colors": {"RNA": ["white", "red"], "sRNA": "magma"},
            "heatmap_interpolation": {"ChIP": "nearest", "mC": "bilinear"}}

    def test_snakefile_translates_both_families(self):
        text = (REPO / "workflow" / "Snakefile").read_text()
        assert 'R_DRAWN_PLOT_FAMILIES = ("heatmaps", "profiles")' in text

    def test_warnings_print_from_the_main_process_only(self):
        # Every cluster job re-parses the Snakefile; the deprecation warning
        # would otherwise land in every job log.
        text = (REPO / "workflow" / "Snakefile").read_text()
        block = text[text.index("PLOT_THEME_JSON, _plot_theme_warnings"):]
        block = block[:block.index("# ----")]
        assert "if workflow.is_main_process:" in block
        assert block.index("if workflow.is_main_process:") < block.index("_plot_theme_warnings:")


def test_colour_map_table_is_sorted_and_complete():
    table = json.loads((REPO / "workflow" / "scripts" / "plot_colormaps.json").read_text())
    maps = table["maps"]
    assert {"seismic", "Oranges", "viridis", "RdBu", "coolwarm"} <= set(maps)
    assert all(len(v) == 17 and all(c.startswith("#") and len(c) == 7 for c in v)
               for v in maps.values())
    assert not any(n.endswith("_r") for n in maps)


def test_profile_rule_hands_the_theme_to_the_r_script():
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

    def test_heatmap_colour_single_name_must_be_a_colour_map(self, tmp_path):
        # One string is a colour-map name; a single colour is not a ramp.
        spec = {"theme": "default", "overrides": {"heatmap_colors": {"RNA": "red"}}}
        proc = _r("", spec, tmp_path)
        assert proc.returncode != 0
        assert "Unknown colour map 'red'" in proc.stderr

    def test_heatmap_colour_by_map_name(self, tmp_path):
        spec = {"theme": "default",
                "overrides": {"heatmap_colors": {"RNA": "viridis", "sRNA": "viridis_r"}}}
        proc = _r("r <- epicc_heatmap_colors(th, 'RNA'); s <- epicc_heatmap_colors(th, 'sRNA'); "
                  "cat(length(r), r[1], r[17], identical(rev(r), s))", spec, tmp_path)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["17", "#440154", "#fde725", "TRUE"]

    def test_colour_map_table_matches_the_seismic_default(self, tmp_path):
        # plot_colormaps.json samples matplotlib; the theme's exact seismic
        # stops are its 0, 1/4, 3/4 and 1 points.
        proc = _r("m <- toupper(epicc_colormap('seismic')); cat(m[c(1, 5, 13, 17)])")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["#00004C", "#0101FF", "#FE0000", "#800000"]

    @pytest.mark.parametrize("dtype,mode", [
        ("ChIP", "bilinear"), ("RNA", "bilinear"), ("mC", "nearest"),
        ("mCHH", "nearest"),        # context -> mC
        ("ChIP_narrow", "bilinear"),
        ("H3K9me2", "bilinear"),    # unknown -> default
    ])
    def test_heatmap_interpolation_defaults_and_fallbacks(self, tmp_path, dtype, mode):
        proc = _r(f"cat(epicc_heatmap_interpolation(th, '{dtype}'))")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.strip() == mode

    def test_heatmap_interpolation_is_an_open_map(self, tmp_path):
        spec = {"theme": "default",
                "overrides": {"heatmap_interpolation": {"mCG": "bilinear", "default": "nearest"}}}
        proc = _r("cat(epicc_heatmap_interpolation(th, 'mCG'), epicc_heatmap_interpolation(th, 'mCHG'), "
                  "epicc_heatmap_interpolation(th, 'ChIP'), epicc_heatmap_interpolation(th, 'H3K4me3'))",
                  spec, tmp_path)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.split() == ["bilinear", "nearest", "bilinear", "nearest"]

    def test_heatmap_interpolation_values(self, tmp_path):
        bad = _r("", {"theme": "default", "overrides": {"heatmap_interpolation": {"RNA": "bicubic"}}},
                 tmp_path)
        assert bad.returncode != 0
        assert "heatmap_interpolation.RNA must be one of bilinear, nearest" in bad.stderr

    def test_heatmap_raster_dpi(self, tmp_path):
        ok = _r("cat(th$heatmap_raster_dpi)",
                {"theme": "default", "overrides": {"heatmap_raster_dpi": 300}}, tmp_path)
        assert ok.returncode == 0, ok.stderr
        assert ok.stdout.strip() == "300"

    def test_profile_type_choices(self, tmp_path):
        ok = _r("cat(th$profile_type)", {"theme": "default", "overrides": {"profile_type": "se"}},
                tmp_path)
        assert ok.returncode == 0, ok.stderr
        assert ok.stdout.strip() == "se"
        bad = _r("", {"theme": "default", "overrides": {"profile_type": "heatmap"}}, tmp_path)
        assert bad.returncode != 0
        assert "profile_type must be one of lines, fill, se, std" in bad.stderr

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
