"""The builder must reject '=' in Levels, as samplefile_validation.py does.

Labels built from levels reach the UpSet scripts as label=type pairs, so '='
is reserved. The builder checks it in three places: the fallback Levels
column (levelsPairIssue, exercised through node), the per-factor level cells,
and factor names, which an imported sheet sets without going through
addFactor.
"""

import os
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest

_REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_BUILDER = os.path.join(_REPO_ROOT, "tools", "epicc-builder.html")
_JS_TEST = os.path.join(os.path.dirname(__file__), "js", "levels_validation_test.js")

sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))
from scripts.samplefile_validation import check_table  # noqa: E402


@pytest.fixture(scope="module")
def html():
    with open(_BUILDER) as fh:
        return fh.read()


class TestBuilderSource:
    def test_add_factor_rejects_equals(self, html):
        m = re.search(r"^function addFactor\(name\)\s*\{(.*?)^\}", html, re.S | re.M)
        assert m, "addFactor not found"
        assert 'name.includes("=")' in m.group(1)

    def test_factor_level_cells_reject_equals(self, html):
        assert 'fval.includes("=")' in html
        assert "level must not contain ':', ',' or '='" in html

    def test_imported_factor_names_are_checked(self, html):
        assert 'factors[fi].includes("=")' in html

    def test_fallback_levels_uses_pair_check(self, html):
        assert "levelsPairIssue(factor, level)" in html


class TestParityWithValidator:
    """Inputs the JS suite rejects are rejected by the pipeline too."""

    @pytest.mark.parametrize("levels", ["genotype:WT=1", "geno=type:WT"])
    def test_validator_rejects(self, levels):
        df = pd.DataFrame([{
            "Sample_ID": "s1", "Assay": "RNAseq", "Genome": "g",
            "Levels": levels, "Replicate_ID": "rep1", "Read_files": "SRR1",
            "Read_layout": "SE", "IP_target": "", "Control": "",
        }])
        with pytest.raises(ValueError, match="must not contain '='"):
            check_table(df)


@pytest.mark.skipif(shutil.which("node") is None, reason="node not available")
def test_js_suite_passes():
    proc = subprocess.run(
        ["node", _JS_TEST], capture_output=True, text=True,
        env={**os.environ, "EPICC_BUILDER_HTML": _BUILDER},
    )
    assert proc.returncode == 0, proc.stdout[-4000:] + "\n" + proc.stderr[-2000:]
    m = re.search(r"(\d+) passed, (\d+) failed", proc.stdout)
    assert m, f"could not read the JS suite summary:\n{proc.stdout[-2000:]}"
    assert int(m.group(1)) >= 7 and int(m.group(2)) == 0
