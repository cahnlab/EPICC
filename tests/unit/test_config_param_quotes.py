"""Quoted tool-parameter strings must be refused at startup.

The builder used to default chip_callpeaks.params to "--keep-dup 'all'
--nomodel"; the quotes survived shell expansion and macs2 rejected the literal
argument 'all'.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), "..", "..", "workflow", "scripts"))

from config_validation import find_quoted_params, format_quoted_params_error


class TestFindQuotedParams:
    def test_flags_the_keep_dup_regression(self):
        cfg = {"chip_callpeaks": {"params": "--keep-dup 'all' --nomodel"}}
        assert find_quoted_params(cfg) == [
            ("chip_callpeaks.params", "--keep-dup 'all' --nomodel")]

    def test_accepts_the_unquoted_form(self):
        cfg = {"chip_callpeaks": {"params": "--keep-dup all --nomodel"}}
        assert find_quoted_params(cfg) == []

    def test_flags_double_quotes_too(self):
        cfg = {"srna_mapping_params": '--mmap "u"'}
        assert find_quoted_params(cfg) == [("srna_mapping_params", '--mmap "u"')]

    @pytest.mark.parametrize("key", [
        "params", "chip_callpeaks_params", "chromap_extra_flags"])
    def test_every_recognised_key_shape(self, key):
        assert find_quoted_params({key: "--x 'y'"}) == [(key, "--x 'y'")]

    def test_nested_under_a_flag_key(self):
        # heatmaps_plot_params maps a label to each flag string.
        cfg = {"heatmaps_plot_params": {"all": "--colorMap 'seismic'"}}
        assert find_quoted_params(cfg) == [
            ("heatmaps_plot_params.all", "--colorMap 'seismic'")]

    def test_ignores_quotes_outside_parameter_settings(self):
        # Only the shell-expanded flag strings are affected; a stray quote in
        # an unrelated setting is not this check's business.
        cfg = {"analysis_name": "Jon's run", "genomes": {"ColCEN": {"genus": "A"}}}
        assert find_quoted_params(cfg) == []

    def test_non_string_flag_values_are_skipped(self):
        # motif_params holds integers.
        assert find_quoted_params({"motif_params": {"n_motifs": 10}}) == []

    def test_reports_every_offender(self):
        cfg = {"chip_callpeaks": {"params": "--a 'b'"},
               "atac_callpeaks": {"params": '--c "d"'}}
        assert len(find_quoted_params(cfg)) == 2


class TestErrorMessage:
    def test_names_the_path_and_the_fix(self):
        msg = format_quoted_params_error(
            [("chip_callpeaks.params", "--keep-dup 'all' --nomodel")])
        assert "chip_callpeaks.params" in msg
        assert "--keep-dup all --nomodel" in msg


class TestShippedOptionsAreClean:
    def test_bundled_options_pass(self):
        yaml = pytest.importorskip("yaml")
        root = os.path.join(os.path.dirname(__file__), "..", "..")
        for name in ("config/epicc-options.yaml",
                     "tests/integration/data/test_options_ChIP.yaml",
                     "tests/integration/data/test_options_ATAC.yaml",
                     "tests/integration/data/test_options_CUT.yaml",
                     "tests/integration/data/test_options_RNA.yaml"):
            path = os.path.join(root, name)
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                cfg = yaml.safe_load(fh)
            assert find_quoted_params(cfg) == [], f"{name} ships quoted params"


class TestBuilderDefaults:
    """The builder must not export a value the pipeline refuses."""

    @staticmethod
    def _param_defaults():
        path = os.path.join(os.path.dirname(__file__), "..", "..",
                            "tools", "epicc-builder.html")
        with open(path) as fh:
            html = fh.read()
        import re
        return re.findall(r'(\w*_?params)\s*:\s*"([^"]*)"', html)

    def test_no_quoted_flag_strings(self):
        offenders = [(k, v) for k, v in self._param_defaults() if "'" in v]
        assert not offenders, f"builder ships quoted params: {offenders}"

    def test_keep_dup_is_unquoted(self):
        values = [v for k, v in self._param_defaults() if "keep-dup" in v]
        assert values, "expected a keep-dup default in the builder"
        for v in values:
            assert "--keep-dup all" in v, v

    def test_defaults_pass_the_startup_guard(self):
        cfg = {k: v for k, v in self._param_defaults()}
        assert find_quoted_params(cfg) == []
