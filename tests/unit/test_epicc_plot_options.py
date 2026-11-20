"""epicc --heatmap-colors / --heatmaps-scales: parsing, and how they reach
the config (heatmaps_scales through --config; heatmap_colors entry by entry
through the merged options file, so the options file's other entries stay)."""

import argparse
import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest
import yaml

_EPICC_SRC = Path(__file__).resolve().parents[2] / "epicc"


def _load_epicc():
    loader = SourceFileLoader("epicc_cli_plot", str(_EPICC_SRC))
    spec = importlib.util.spec_from_loader("epicc_cli_plot", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli_plot", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


def _parse(*argv):
    p = argparse.ArgumentParser()
    epicc.add_plot_options(p)
    return p.parse_args(list(argv))


class TestParse:
    def test_maps_and_colour_lists(self):
        got = epicc._parse_heatmap_colors("H3K9me2=Reds, RNA=white:darkgreen,RAMPAGE=PuOr_r")
        assert got == {"H3K9me2": "Reds", "RNA": ["white", "darkgreen"], "RAMPAGE": "PuOr_r"}

    def test_hex_colours(self):
        assert epicc._parse_heatmap_colors("mCG=#ffffff:#ff0000") == \
            {"mCG": ["#ffffff", "#ff0000"]}

    @pytest.mark.parametrize("bad", ["H3K9me2", "=Reds", "H3K9me2=", "RNA=white::red", ","])
    def test_malformed(self, bad):
        with pytest.raises(argparse.ArgumentTypeError):
            epicc._parse_heatmap_colors(bad)

    def test_scales_choices(self):
        assert _parse("--heatmaps-scales", "shared").heatmaps_scales == "shared"
        with pytest.raises(SystemExit):
            _parse("--heatmaps-scales", "global")

    def test_unset_by_default(self):
        args = _parse()
        assert args.heatmap_colors is None and args.heatmaps_scales is None


class TestReachesTheConfig:
    def test_scales_go_through_config(self):
        args = _parse("--heatmaps-scales", "shared")
        assert "heatmaps_scales=shared" in epicc._collect_config_overrides(args)

    def test_colours_merge_entry_by_entry(self, tmp_path, monkeypatch):
        opts = tmp_path / "opts.yaml"
        opts.write_text(yaml.safe_dump({
            "heatmaps_scales": "type",
            "plot_theme_overrides": {"base_size": 12,
                                     "heatmap_colors": {"mC": "Greens", "RNA": "viridis"}}}))
        monkeypatch.chdir(tmp_path)
        args = _parse("--heatmap-colors", "H3K9me2=Reds,RNA=white:darkgreen")
        args.options = str(opts)
        args.te_file = None
        merged = yaml.safe_load(Path(epicc._prepare_options_file(args)).read_text())
        assert merged["plot_theme_overrides"] == {
            "base_size": 12,
            "heatmap_colors": {"mC": "Greens", "RNA": ["white", "darkgreen"], "H3K9me2": "Reds"}}

    def test_colours_without_an_overrides_block(self, tmp_path, monkeypatch):
        opts = tmp_path / "opts.yaml"
        opts.write_text("plot_theme: default\nplot_theme_overrides:\n")
        monkeypatch.chdir(tmp_path)
        args = _parse("--heatmap-colors", "ChIP=Blues")
        args.options = str(opts)
        args.te_file = None
        merged = yaml.safe_load(Path(epicc._prepare_options_file(args)).read_text())
        assert merged["plot_theme_overrides"] == {"heatmap_colors": {"ChIP": "Blues"}}

    @pytest.mark.parametrize("parser", ["p_run", "p_val", "p_out"])
    def test_run_validate_and_output_take_them(self, parser):
        assert f"add_plot_options({parser})" in _EPICC_SRC.read_text()
