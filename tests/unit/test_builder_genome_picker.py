"""Tests for the builder's multi-select Genome column editor.

The Genome cell holds a comma-separated list of references ("B73,W22"), which a
single-entry dropdown cannot express, so the column uses ``genomeMultiEditor``.

Source checks run off the HTML alone. The behavioral suite runs through node
(``tests/unit/js/genome_picker_test.js``) and is skipped where node is absent.
"""

import os
import re
import shutil
import subprocess
import sys

import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
_BUILDER = os.path.join(_REPO_ROOT, "tools", "epicc-builder.html")
_JS_TEST = os.path.join(os.path.dirname(__file__), "js", "genome_picker_test.js")

sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))

from sample_sheet import read_sample_sheet  # noqa: E402


@pytest.fixture(scope="module")
def html():
    with open(_BUILDER) as fh:
        return fh.read()


class TestBuilderSource:
    """Checks that need only the HTML source."""

    def _genome_column(self, html):
        m = re.search(
            r'title: "Genome", field: "Genome".*?\}, descHooks\("Genome"\)\)',
            html, re.S,
        )
        assert m, "Genome column definition not found in the builder"
        return m.group(0)

    def test_genome_column_uses_the_multiselect_editor(self, html):
        assert "editor: genomeMultiEditor" in self._genome_column(html)

    def test_genome_column_is_not_a_single_entry_list(self, html):
        """Tabulator's list editor holds one value, making multi-genome rows
        unenterable through the UI even though the loader still accepts them.
        """
        col = self._genome_column(html)
        assert 'editor: "list"' not in col
        assert "freetext" not in col

    def test_editor_and_helper_are_defined(self, html):
        assert re.search(r"^function genomeMultiEditor\(", html, re.M)
        assert re.search(r"^function parseGenomeList\(", html, re.M)

    def test_picker_styles_exist(self, html):
        """The panel is appended to document.body, so it needs its own styles;
        without them it renders unstyled at the top-left of the page.
        """
        for cls in (
            ".genome-picker",
            ".genome-picker-item",
            ".genome-picker-add",
            ".genome-picker-empty",
        ):
            assert cls + " {" in html, f"missing CSS for {cls}"

    def test_picker_sits_above_the_context_menu(self, html):
        """Right-clicking a row opens .ctx-menu at z-index 10000."""
        m = re.search(r"\.genome-picker \{(.*?)\}", html, re.S)
        assert m
        z = re.search(r"z-index:\s*(\d+)", m.group(1))
        assert z and int(z.group(1)) > 10000, "picker must stack above .ctx-menu"

    def test_column_description_mentions_multiple_genomes(self, html):
        """The hover description is the only in-app explanation of the feature."""
        m = re.search(r"^\s*Genome: \"(.*?)\",$", html, re.M | re.S)
        assert m, "Genome column description not found"
        desc = m.group(1)
        assert "Tick" in desc, "description should explain the tick-box picker"
        assert "," in desc and "B73" in desc, "should show a comma-separated example"


class TestParseParityWithLoader:
    """The builder and the loader must agree on how a Genome cell is parsed.

    If they drift, the builder shows a sheet as valid that the pipeline reads
    differently (or vice versa).
    """

    HEADER = ("Assay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
              "Read_layout\tIP_target\tControl")

    def _genomes_for(self, tmp_path, cell):
        f = tmp_path / "s.tsv"
        f.write_text(
            "\n".join([
                self.HEADER,
                f"ATAC\t{cell}\tgenotype:WT\trep1\tSRR1\tSE\t\t",
            ]) + "\n"
        )
        return list(read_sample_sheet(f)["Genome"])

    @pytest.mark.parametrize("cell,expected", [
        ("B73,W22", ["B73", "W22"]),
        (" B73 , W22 ", ["B73", "W22"]),
        ("B73,,W22", ["B73", "W22"]),      # stray comma tolerated on both sides
        ("ColCEN", ["ColCEN"]),
    ])
    def test_loader_matches_the_js_cases(self, tmp_path, cell, expected):
        """Same inputs the JS suite asserts parseGenomeList on."""
        assert self._genomes_for(tmp_path, cell) == expected


@pytest.mark.skipif(shutil.which("node") is None, reason="node not available")
class TestPickerBehavior:
    def test_js_suite_passes(self):
        proc = subprocess.run(
            ["node", _JS_TEST],
            capture_output=True, text=True,
            env={**os.environ, "EPICC_BUILDER_HTML": os.path.abspath(_BUILDER)},
        )
        assert proc.returncode == 0, (
            "genome picker JS suite failed:\n"
            + proc.stdout[-4000:] + "\n" + proc.stderr[-2000:]
        )
        # Guard against the suite silently becoming a no-op (e.g. a rename makes
        # extract() return nothing and every assertion vanishes).
        m = re.search(r"(\d+) passed, (\d+) failed", proc.stdout)
        assert m, f"could not read the JS suite summary:\n{proc.stdout[-2000:]}"
        assert int(m.group(1)) >= 25, f"only {m.group(1)} JS assertions ran"
        assert int(m.group(2)) == 0
