"""Builder-source tests for issue #66.

Two independent complaints, one shared root cause for the ordering half:

1. A ChIP Input used only as a control was rejected for having no Peak_type.
   The builder's rule must match `samplefile_validation.py`: required only on
   a pulldown row that declares a Control.

2. Row order was not preserved through export or import. `descHooks` attached
   `headerClick` to every data column to show its description, and Tabulator
   runs that callback *alongside* its built-in header sort — so clicking a
   header to read the docs silently re-sorted the table. `table.getData()`
   then handed that order to the exported sheet, and an active sort survives
   `setData`, so an import was re-sorted too.

These are source-level checks: the code under test is welded to a live
Tabulator instance, so behavioural tests would need a browser. The parity test
against the Python validator is the one that would catch real drift.
"""

import os
import re
import sys

import pandas as pd
import pytest

_REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_BUILDER = os.path.join(_REPO_ROOT, "tools", "epicc-builder.html")

sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))
from scripts.samplefile_validation import check_table  # noqa: E402


@pytest.fixture(scope="module")
def html():
    with open(_BUILDER) as fh:
        return fh.read()


def _desc_hooks_body(html):
    m = re.search(r"function descHooks\(field\)\s*\{(.*?)\n\}", html, re.S)
    assert m, "descHooks() not found in the builder"
    return m.group(1)


class TestHeaderClickNoLongerSorts:
    def test_desc_hooks_has_no_header_click(self, html):
        assert "headerClick" not in _desc_hooks_body(html)

    def test_hover_and_cell_click_still_show_descriptions(self, html):
        body = _desc_hooks_body(html)
        assert "headerMouseEnter" in body
        assert "cellClick" in body

    def test_no_column_reattaches_header_click(self, html):
        # A single stray headerClick on any column reintroduces the accidental
        # sort for that column -- including the dynamic factor columns, which
        # build their hooks inline rather than through descHooks(). Comment
        # lines are stripped so the prose explaining this rule doesn't match.
        code = [l for l in html.split("\n") if not l.lstrip().startswith("//")]
        offenders = [l.strip() for l in code if "headerClick" in l]
        assert offenders == [], offenders


class TestImportPreservesFileOrder:
    def test_import_clears_any_active_sort(self, html):
        m = re.search(r"function importTSV\(text\)\s*\{(.*?)\n\}", html, re.S)
        assert m, "importTSV() not found"
        body = m.group(1)
        assert "clearSort" in body, "importTSV must clear the sort before loading"
        assert body.index("clearSort") < body.index("table.setData"), \
            "the sort must be cleared before setData, or it is re-applied"

    def test_export_writes_rows_in_table_order(self, html):
        # buildExportMatrix must walk table.getData() straight through -- no
        # sort, reverse, or grouping between the table and the sheet.
        m = re.search(r"function buildExportMatrix\(\)\s*\{(.*?)\n\}", html, re.S)
        assert m, "buildExportMatrix() not found"
        body = m.group(1)
        for forbidden in (".sort(", ".reverse(", "groupBy"):
            assert forbidden not in body, f"export must not {forbidden} rows"


class TestPeakTypeRuleMatchesValidator:
    """The builder's Peak_type rule must agree with the Python validator.

    Drift here is the actual hazard: a builder that accepts a sheet the
    pipeline rejects (or vice versa) sends the user round a loop.
    """

    def test_builder_requires_peak_type_only_with_a_control(self, html):
        m = re.search(r"if \(PEAK_TYPE_ASSAYS\.has\(assay\)\) \{(.*?)\n    \} else if", html, re.S)
        assert m, "the Peak_type validation branch was not found"
        body = m.group(1)
        assert "!peakType && control" in body, \
            "blank Peak_type must only be an error when a Control is set"

    @pytest.mark.parametrize("ip_pt,ctrl_pt,should_pass", [
        ("broad", "",       True),   # control-only row may omit it
        ("broad", "narrow", True),   # ... and may still carry one
        ("",      "broad",  False),  # a peak caller may not omit it
        ("wide",  "",       False),  # invalid value on the caller
        ("broad", "wide",   False),  # invalid value on the control row
    ])
    def test_validator_agrees_with_that_rule(self, tmp_path, ip_pt, ctrl_pt, should_pass):
        ip_fq = tmp_path / "ip.fq.gz"
        ip_fq.write_bytes(b"")
        ctrl_fq = tmp_path / "ctrl.fq.gz"
        ctrl_fq.write_bytes(b"")

        def row(sid, ip_target, control, peak_type, reads):
            return {
                "Sample_ID": sid, "Assay": "ChIP", "Genome": "test_genome",
                "Levels": "genotype:WT", "Replicate_ID": "rep1",
                "Read_files": str(reads), "Read_layout": "SE",
                "IP_target": ip_target, "Control": control,
                "Peak_type": peak_type,
            }

        df = pd.DataFrame([
            row("ip1", "H3K9me2", "ctrl1", ip_pt, ip_fq),
            row("ctrl1", "Input", "", ctrl_pt, ctrl_fq),
        ])
        if should_pass:
            check_table(df)
        else:
            with pytest.raises(ValueError):
                check_table(df)

    def test_description_documents_the_relaxed_rule(self, html):
        m = re.search(r"Peak_type: \"Peak-calling type[^\"]*\"", html)
        assert m, "Peak_type column description not found"
        assert "Control" in m.group(0), \
            "the description should say the requirement follows Control"


class TestShippedExamplesFollowTheRule:
    def test_builder_example_control_rows_have_no_peak_type(self, html):
        # The shipped example is how most users learn the format.
        for sid in ("WT_leaf_Input_rep1", "WT_leaf_Input_rep2"):
            m = re.search(r'Sample_ID: "%s", Assay: "ChIP", Peak_type: "([^"]*)"' % sid, html)
            assert m, f"example row {sid} not found"
            assert m.group(1) == "", f"{sid} is control-only; Peak_type should be blank"

    def test_example_samples_tsv_control_rows_have_no_peak_type(self):
        path = os.path.join(_REPO_ROOT, "config", "example_samples.tsv")
        with open(path) as fh:
            lines = [l.rstrip("\n") for l in fh]
        hi = next(i for i, l in enumerate(lines) if l.startswith("Sample_ID\t"))
        cols = {h: i for i, h in enumerate(lines[hi].split("\t"))}
        checked = 0
        for line in lines[hi + 1:]:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) <= cols["Peak_type"]:
                continue
            if f[cols["Assay"]] not in ("ChIP", "CUT_RUN", "CUT_TAG"):
                continue
            if f[cols["Control"]].strip():
                continue
            assert f[cols["Peak_type"]].strip() == "", (
                f"{f[cols['Sample_ID']]} is control-only; Peak_type should be blank"
            )
            checked += 1
        assert checked > 0, "no control-only pulldown rows found to check"
