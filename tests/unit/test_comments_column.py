"""Unit tests for the free-text Comments column (issue #52).

Comments is a user annotation the pipeline never reads: optional in both
directions (may be absent from a sheet, may be blank on a row), unconstrained,
and inert — it must never leak into analysis keys, names, merge keys, or the
builder's assay-gating.
"""

import os
import re
import sys

import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))

from sample_sheet import (  # noqa: E402
    NEW_COLNAMES,
    OPTIONAL_COLNAMES,
    SHEET_COLNAMES,
    add_compat_columns,
    build_analysis_key,
    build_analysis_name,
    build_control_merge_key,
    read_sample_sheet,
)

_FULL_HEADER = (
    "Assay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
    "Read_layout\tIP_target\tControl\tPeak_type\tComments"
)
# The pre-Comments header, for the back-compat case.
_LEGACY_HEADER = (
    "Assay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
    "Read_layout\tIP_target\tControl\tPeak_type"
)


def _sheet(tmp_path, header, rows, name="s.tsv"):
    f = tmp_path / name
    f.write_text("\n".join([header] + rows) + "\n")
    return read_sample_sheet(f)


class TestSchemaPosition:
    def test_comments_is_last_column(self):
        assert NEW_COLNAMES[-1] == "Comments"
        assert NEW_COLNAMES[-2] == "Peak_type"

    def test_comments_is_optional(self):
        assert "Comments" in OPTIONAL_COLNAMES

    def test_builder_columns_match_sheet_colnames(self):
        """The builder's COLUMNS and SHEET_COLNAMES must not drift apart.

        Both express the canonical exported column order; #55 had to keep five
        places in step and this is the one pair a test can hold together
        directly. Neither holds Sample_ID: it is derived, never exported (#91).
        """
        html = open(
            os.path.join(_REPO_ROOT, "tools", "epicc-builder.html")
        ).read()
        m = re.search(r"const COLUMNS = \[(.*?)\];", html, re.S)
        assert m, "COLUMNS array not found in the builder"
        cols = re.findall(r'"([^"]+)"', m.group(1))
        assert cols == SHEET_COLNAMES, (
            "builder COLUMNS and SHEET_COLNAMES disagree:\n"
            f"  builder: {cols}\n  python : {SHEET_COLNAMES}"
        )


class TestOptionalAndRoundTrip:
    ROW = ("RNAseq\tG\tgenotype:WT\trep1\tSRR1\tSE\t\t\t")

    def test_absent_column_still_loads(self, tmp_path):
        # A sheet written before the column existed must parse, with "" filled in.
        df = _sheet(tmp_path, _LEGACY_HEADER,
                    ["RNAseq\tG\tgenotype:WT\trep1\tSRR1\tSE\t\t\t"])
        assert "Comments" in df.columns
        assert df["Comments"].tolist() == [""]

    def test_blank_comment_is_empty_string_not_nan(self, tmp_path):
        df = _sheet(tmp_path, _FULL_HEADER, [self.ROW + "\t"])
        assert df["Comments"].tolist() == [""]

    def test_freetext_with_spaces_and_punctuation(self, tmp_path):
        note = "batch 3; re-prepped 2026-04-01 (low yield), see notebook p.42!"
        df = _sheet(tmp_path, _FULL_HEADER, [self.ROW + "\t" + note])
        assert df["Comments"].iloc[0] == note

    def test_hash_within_field_is_preserved(self, tmp_path):
        """A '#' mid-field must survive: read_sample_sheet drops only
        FULL-LINE comments, deliberately not using pandas' comment='#'."""
        note = "lane #3 failed QC, #rerun"
        df = _sheet(tmp_path, _FULL_HEADER, [self.ROW + "\t" + note])
        assert df["Comments"].iloc[0] == note

    def test_full_line_comment_rows_still_skipped(self, tmp_path):
        df = _sheet(tmp_path, _FULL_HEADER, [
            "# parked sample, do not run",
            self.ROW + "\tkeep me",
        ])
        assert df["Sample_ID"].tolist() == ["WT_RNA_rep1"]
        assert df["Comments"].iloc[0] == "keep me"


class TestInertness:
    """Comments must not influence any derived identity."""

    BASE = "ChIP\tG\tgenotype:WT\trep1\tSRR1\tSE\tH3K9me2\tWT_Input_rep1\tbroad"

    def _row(self, tmp_path, comment, name):
        df = _sheet(tmp_path, _FULL_HEADER, [self.BASE + "\t" + comment], name)
        return add_compat_columns(df).iloc[0]

    def test_analysis_name_unaffected(self, tmp_path):
        a = self._row(tmp_path, "", "a.tsv")
        b = self._row(tmp_path, "some note about the batch", "b.tsv")
        assert build_analysis_name(a) == build_analysis_name(b)
        assert build_analysis_name(a) == "ChIP_broad__WT__H3K9me2__G"

    def test_sample_id_unaffected(self, tmp_path):
        # Sample_ID is derived from the row (#91); Comments must not feed it.
        a = self._row(tmp_path, "", "a.tsv")
        b = self._row(tmp_path, "some note about the batch", "b.tsv")
        assert a["Sample_ID"] == b["Sample_ID"] == "WT_H3K9me2_rep1"

    def test_analysis_key_unaffected(self, tmp_path):
        a = self._row(tmp_path, "", "a.tsv")
        b = self._row(tmp_path, "note", "b.tsv")
        assert build_analysis_key(a) == build_analysis_key(b)

    def test_control_merge_key_unaffected(self, tmp_path):
        a = self._row(tmp_path, "", "a.tsv")
        b = self._row(tmp_path, "note", "b.tsv")
        assert build_control_merge_key(a) == build_control_merge_key(b)

    def test_not_read_anywhere_in_workflow_logic(self):
        """Guard against someone wiring Comments into pipeline behaviour.

        Only the schema definitions in sample_sheet.py may mention it.
        """
        offenders = []
        for sub in ("rules", "scripts"):
            base = os.path.join(_REPO_ROOT, "workflow", sub)
            for dirpath, _, filenames in os.walk(base):
                for fn in filenames:
                    if not fn.endswith((".py", ".smk")):
                        continue
                    path = os.path.join(dirpath, fn)
                    for i, line in enumerate(open(path), 1):
                        if "Comments" not in line:
                            continue
                        # schema declaration + its explanatory comments are fine
                        if fn == "sample_sheet.py":
                            continue
                        offenders.append(f"{fn}:{i}: {line.strip()}")
        assert not offenders, (
            "Comments must stay inert; found references outside the schema:\n"
            + "\n".join(offenders)
        )


class TestBuilderWiring:
    """Checks against the builder source that don't need a browser."""

    @pytest.fixture(scope="class")
    def html(self):
        return open(os.path.join(_REPO_ROOT, "tools", "epicc-builder.html")).read()

    def test_comments_not_assay_gated(self, html):
        m = re.search(r"const ASSAY_GATED_FIELDS = \{(.*?)\};", html, re.S)
        assert m, "ASSAY_GATED_FIELDS not found"
        assert "Comments" not in m.group(1), (
            "Comments must be editable on every row, not assay-gated"
        )

    def test_comments_optional_on_import(self, html):
        m = re.search(r"var OPTIONAL_ON_IMPORT = \[(.*?)\];", html, re.S)
        assert m, "OPTIONAL_ON_IMPORT not found"
        assert "Comments" in m.group(1)

    def test_sanitize_strips_tabs_and_newlines(self, html):
        # The export sanitizer is what keeps free text from breaking the TSV.
        m = re.search(r"function sanitizeCell\(v\) \{(.*?)\n\}", html, re.S)
        assert m, "sanitizeCell not found"
        assert re.search(r"replace\(/\[\\t\\r\\n\]\+/g", m.group(1)), (
            "sanitizeCell no longer strips tab/CR/LF — free-text Comments "
            "could break the exported TSV"
        )
