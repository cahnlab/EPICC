"""Tests for commented samples in the builder (issue #65).

A sample-sheet line whose first non-whitespace character is '#' is skipped by
``read_sample_sheet``. The builder imports such a line as a greyed-out
"commented" row, and writes it back with its '#', so a round trip is lossless.

The behavioral suite runs through node (``tests/unit/js/commented_rows_test.js``)
and is skipped where node is absent.
"""

import json
import os
import re
import shutil
import subprocess
import sys

import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
_BUILDER = os.path.join(_REPO_ROOT, "tools", "epicc-builder.html")
_JS_TEST = os.path.join(os.path.dirname(__file__), "js", "commented_rows_test.js")

sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))

from sample_sheet import read_sample_sheet  # noqa: E402

_needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node not available")

_HEADER = (
    "Sample_ID\tAssay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
    "Read_layout\tIP_target\tControl\tPeak_type\tComments"
)


def _row(sid, rep, srr, prefix=""):
    return prefix + "\t".join([
        sid, "ATAC", "ColCEN", "genotype:WT", rep, srr, "SE", "", "", "", "",
    ])


@pytest.fixture(scope="module")
def html():
    with open(_BUILDER) as fh:
        return fh.read()


class TestBuilderSource:
    def test_toolbar_has_comment_action_next_to_duplicate_and_delete(self, html):
        toolbar = re.search(r'id="btn-dup-selected".*?id="btn-del-selected"', html, re.S)
        assert toolbar and 'id="btn-comment-selected"' in toolbar.group(0)

    def test_commented_rows_are_styled(self, html):
        assert ".row-commented" in html
        assert 'classList.toggle("row-commented"' in html


@_needs_node
class TestBehavior:
    def test_node_suite(self):
        proc = subprocess.run(["node", _JS_TEST], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stdout + proc.stderr

    def test_builder_and_pipeline_agree_on_active_samples(self, tmp_path):
        """Every row the builder treats as active is one the pipeline reads, and
        every row it flags commented is one the pipeline skips."""
        sheet = tmp_path / "s.tsv"
        sheet.write_text("\n".join([
            "# preamble",
            _HEADER,
            _row("a_rep1", "rep1", "SRR1"),
            "",
            "# a note",
            _row("b_rep1", "rep1", "SRR2", prefix="#"),
            _row("c_rep1", "rep1", "SRR3", prefix="  # "),
            "#" + _HEADER,
            _row("d_rep1", "rep1", "SRR4"),
        ]) + "\n")
        proc = subprocess.run(["node", _JS_TEST, str(sheet)], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        view = json.loads(proc.stdout)
        active = sorted(r["Sample_ID"] for r in view if not r["commented"])
        parked = sorted(r["Sample_ID"] for r in view if r["commented"])
        assert parked == ["b_rep1", "c_rep1"]
        assert active == sorted(read_sample_sheet(sheet)["Sample_ID"].unique())
