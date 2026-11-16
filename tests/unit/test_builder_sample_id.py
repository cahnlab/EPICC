"""Sample_ID is derived from the other columns, never read from the sheet (#91).

Covers derive_sample_id, the rejection of a sheet that still has a Sample_ID
column, the validation messages that replace the old Sample_ID checks, and
that the builder (tools/epicc-builder.html, run through node) derives the same
names as the pipeline for every shipped sheet.
"""

import glob
import json
import os
import shutil
import subprocess
import sys

import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))

from sample_sheet import (  # noqa: E402
    SHEET_COLNAMES,
    add_compat_columns,
    derive_sample_id,
    read_sample_sheet,
)
from scripts.samplefile_validation import check_table  # noqa: E402

_JS_VIEW = os.path.join(os.path.dirname(__file__), "js", "commented_rows_test.js")
_needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node not available")

HEADER = "\t".join(SHEET_COLNAMES)


def _r(assay, levels, rep, ip="", control="", peak="", genome="G", reads="SRR1",
       layout="SE", comments=""):
    return "\t".join([assay, genome, levels, rep, reads, layout, ip, control,
                      peak, comments])


def _read(tmp_path, rows, header=HEADER, name="s.tsv"):
    f = tmp_path / name
    f.write_text("\n".join([header] + rows) + "\n")
    return read_sample_sheet(f)


class TestDerivation:
    @pytest.mark.parametrize("row,expected", [
        ({"Assay": "ChIP", "Levels": "genotype:WT", "IP_target": "H3K9me2",
          "Replicate_ID": "rep1"}, "WT_H3K9me2_rep1"),
        ({"Assay": "ChIP", "Levels": "genotype:WT,tissue:leaf", "IP_target": "Input",
          "Replicate_ID": "rep2"}, "WT_leaf_Input_rep2"),
        ({"Assay": "CUT_RUN", "Levels": "genotype:WT", "IP_target": "H3K27me3",
          "Replicate_ID": "rep1"}, "WT_CR_H3K27me3_rep1"),
        ({"Assay": "CUT_TAG", "Levels": "genotype:WT", "IP_target": "IgG",
          "Replicate_ID": "rep1"}, "WT_CT_IgG_rep1"),
        ({"Assay": "RNAseq", "Levels": "genotype:WT", "Replicate_ID": "rep1"}, "WT_RNA_rep1"),
        ({"Assay": "WGBS_nd", "Levels": "genotype:WT", "Replicate_ID": "rep1"}, "WT_WGBS_nd_rep1"),
        ({"Assay": "dmC", "Levels": "tissue:cell_line", "Replicate_ID": "HG002"},
         "cell_line_dmC_HG002"),
        # The combined internal token names the library the same way.
        ({"Assay": "ChIP_narrow", "Levels": "genotype:WT", "IP_target": "CTCF",
          "Replicate_ID": "rep1"}, "WT_CTCF_rep1"),
    ])
    def test_forms(self, row, expected):
        assert derive_sample_id(row) == expected

    @pytest.mark.parametrize("row", [
        {"Assay": "ChIP", "Levels": "genotype:WT", "Replicate_ID": "rep1"},  # no IP
        {"Assay": "RNAseq", "Levels": "genotype:WT", "Replicate_ID": ""},
        {"Assay": "RNAseq", "Levels": "", "Replicate_ID": "rep1"},
        {"Assay": "RNAseq", "Levels": "WT", "Replicate_ID": "rep1"},  # malformed
        {"Assay": "Nope", "Levels": "genotype:WT", "Replicate_ID": "rep1"},
    ])
    def test_underivable_is_blank(self, row):
        assert derive_sample_id(row) == ""

    @pytest.mark.parametrize("field,value", [
        ("Genome", "Other"), ("Peak_type", "narrow"), ("Read_files", "SRR9"),
        ("Read_layout", "PE"), ("Control", "x"), ("Comments", "note"),
    ])
    def test_non_identity_fields_do_not_contribute(self, field, value):
        base = {"Assay": "ChIP", "Levels": "genotype:WT", "IP_target": "H3",
                "Replicate_ID": "rep1", "Genome": "G", "Peak_type": "broad",
                "Read_files": "SRR1", "Read_layout": "SE", "Control": "",
                "Comments": ""}
        assert derive_sample_id({**base, field: value}) == derive_sample_id(base)


class TestReadSampleSheet:
    def test_no_sample_id_column_needed(self, tmp_path, capsys):
        df = _read(tmp_path, [_r("RNAseq", "genotype:WT", "rep1")])
        assert df["Sample_ID"].tolist() == ["WT_RNA_rep1"]
        assert "Sample_ID column" not in capsys.readouterr().out

    def test_multi_genome_row_shares_one_id(self, tmp_path):
        df = add_compat_columns(_read(tmp_path, [_r("RNAseq", "genotype:WT", "rep1",
                                                    genome="B73,W22")]))
        assert set(df["Sample_ID"]) == {"WT_RNA_rep1"}
        assert sorted(df["mapped_name"]) == ["WT_RNA_rep1__B73", "WT_RNA_rep1__W22"]

    def test_sample_id_column_is_rejected(self, tmp_path):
        # No compatibility layer: a sheet that still names its samples is an
        # error, not silently reinterpreted.
        with pytest.raises(ValueError, match="has a Sample_ID column"):
            _read(tmp_path, ["WT_RNA_rep1\t" + _r("RNAseq", "genotype:WT", "rep1")],
                  header="Sample_ID\t" + HEADER)


class TestValidation:
    def test_duplicate_derived_id_explains_the_fix(self, tmp_path):
        df = _read(tmp_path, [
            _r("RNAseq", "genotype:WT", "rep1", reads="SRR1"),
            _r("RNAseq", "genotype:WT", "rep1", reads="SRR2"),
        ])
        with pytest.raises(ValueError) as e:
            check_table(df, check_paths=False)
        msg = str(e.value)
        assert "same Sample_ID as row #1" in msg and "joined with '+'" in msg

    def test_same_id_on_different_genomes_is_allowed(self, tmp_path):
        df = _read(tmp_path, [
            _r("RNAseq", "genotype:WT", "rep1", genome="B73"),
            _r("RNAseq", "genotype:WT", "rep1", genome="W22"),
        ])
        check_table(df, check_paths=False)

    def test_underivable_row_reported(self, tmp_path):
        df = _read(tmp_path, [_r("ChIP", "genotype:WT", "rep1", ip="")])
        with pytest.raises(ValueError, match="Sample_ID could not be derived"):
            check_table(df, check_paths=False)

    @pytest.mark.parametrize("kw,col", [
        ({"rep": "rep 1"}, "Replicate_ID"),
        ({"ip": "H3/K9"}, "IP_target"),
        ({"levels": "genotype:W T"}, "Levels"),
    ])
    def test_unsafe_component_names_its_column(self, tmp_path, kw, col):
        args = {"assay": "ChIP", "levels": "genotype:WT", "rep": "rep1", "ip": "H3"}
        args.update(kw)
        df = _read(tmp_path, [_r(args["assay"], args["levels"], args["rep"], ip=args["ip"])])
        with pytest.raises(ValueError, match=f"{col} value .* unsafe characters"):
            check_table(df, check_paths=False)

    @pytest.mark.parametrize("rep", ["_rep1", "rep1_", "rep__1"])
    def test_component_that_would_make_double_underscore(self, tmp_path, rep):
        df = _read(tmp_path, [_r("RNAseq", "genotype:WT", rep)])
        with pytest.raises(ValueError, match="Replicate_ID value .* must not contain '__'"):
            check_table(df, check_paths=False)

    def test_unknown_control_suggests_close_match(self, tmp_path):
        df = _read(tmp_path, [
            _r("ChIP", "genotype:WT", "rep1", ip="Input", reads="SRR1"),
            _r("ChIP", "genotype:WT", "rep1", ip="H3", control="WT_Input_rep2x",
               peak="broad", reads="SRR2"),
        ])
        with pytest.raises(ValueError, match="did you mean 'WT_Input_rep1'"):
            check_table(df, check_paths=False)


_SHEETS = sorted(
    [os.path.join(_REPO_ROOT, "config", "example_samples.tsv")]
    + glob.glob(os.path.join(_REPO_ROOT, "tests", "integration", "data", "test_samples_*.tsv"))
)


@pytest.mark.parametrize("sheet", _SHEETS, ids=os.path.basename)
class TestShippedSheets:
    def test_no_sample_id_column(self, sheet):
        with open(sheet) as fh:
            header = next(l for l in fh if l.strip() and not l.lstrip().startswith("#"))
        assert "Sample_ID" not in header.split("\t")

    def test_validates(self, sheet):
        check_table(read_sample_sheet(sheet), check_paths=False)

    @_needs_node
    def test_builder_derives_the_pipelines_ids(self, sheet):
        proc = subprocess.run(["node", _JS_VIEW, sheet], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        builder = sorted(r["Sample_ID"] for r in json.loads(proc.stdout)
                         if not r["commented"])
        assert builder == sorted(read_sample_sheet(sheet)["Sample_ID"].unique())
