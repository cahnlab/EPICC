"""Per-replicate count tables must be addressed by mapped_name.

The counts are a post-alignment product, so they carry the genome. Both rules
here are `run:` blocks, which a dry-run never executes, so a mismatch only
surfaces during a real run (#71 for RNA, #79 for sRNA).
"""

import os
import re
import sys

import pytest

_REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_RNASEQ_SMK = os.path.join(_REPO_ROOT, "workflow", "rules", "RNAseq.smk")
_SRNA_SMK = os.path.join(_REPO_ROOT, "workflow", "rules", "smallRNA.smk")

sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))
from scripts.sample_sheet import (  # noqa: E402
    add_compat_columns, read_sample_sheet,
)


@pytest.fixture(scope="module")
def srna_prep_rule():
    with open(_SRNA_SMK) as fh:
        src = fh.read()
    m = re.search(r"rule prep_files_for_differential_srna_clusters:(.*?)(?=\nrule )",
                  src, re.S)
    assert m, "prep_files_for_differential_srna_clusters not found"
    return m.group(1)


@pytest.fixture(scope="module")
def srna_grouped_input():
    with open(_SRNA_SMK) as fh:
        src = fh.read()
    m = re.search(r"def define_input_for_grouped_analysis\(.*?\n(?=\S)", src, re.S)
    assert m, "define_input_for_grouped_analysis not found"
    return m.group(0)


@pytest.fixture(scope="module")
def deg_rule():
    with open(_RNASEQ_SMK) as fh:
        src = fh.read()
    m = re.search(r"rule prep_files_for_DEGs:(.*?)(?=\nrule )", src, re.S)
    assert m, "prep_files_for_DEGs not found"
    return m.group(1)


class TestCountsPathUsesMappedName:
    def test_counts_path_is_built_from_mapped_name(self, deg_rule):
        m = re.search(r'counts__\{(\w+)\}\.tab', deg_rule)
        assert m, "the per-replicate counts path was not found"
        assert m.group(1) == "mname", (
            f"counts path uses {{{m.group(1)}}}; it must use the mapped_name")

    def test_replicate_frame_selects_mapped_name(self, deg_rule):
        assert "filtered_samples[['mapped_name', 'Replicate']]" in deg_rule, \
            "the replicate frame must carry mapped_name, not sample_name"

    def test_bare_sample_name_is_not_used_for_counts(self, deg_rule):
        assert "counts__{sname}.tab" not in deg_rule


class TestNamesActuallyDiffer:
    """Guards the assumption the fix rests on: the two names are not equal."""

    def test_mapped_name_carries_the_genome(self):
        sheet = os.path.join(_REPO_ROOT, "tests", "integration", "data",
                             "test_samples_pombe.tsv")
        df = add_compat_columns(read_sample_sheet(sheet))
        rna = df[df["Assay"] == "RNAseq"]
        assert not rna.empty, "no RNAseq rows in the pombe fixture"
        for _, row in rna.iterrows():
            assert row["mapped_name"] != row["sample_name"]
            assert row["mapped_name"].endswith("__" + row["Genome"])


class TestShortStackCountsColumns:
    """ShortStack names its count columns after the BAM basenames it was given.

    define_input_for_grouped_analysis builds those from mapped_name, so the
    column selection has to use the same name (#79).
    """

    def test_bam_inputs_are_built_from_mapped_name(self, srna_grouped_input):
        assert "clean__{mname}_condensed.bam" in srna_grouped_input, (
            "the ShortStack BAM path must be built from mapped_name")

    def test_column_order_uses_mapped_name(self, srna_prep_rule):
        assert "ROW['mapped_name']" in srna_prep_rule, \
            "count columns must be selected by mapped_name"
        assert "ROW['sample_name']" not in srna_prep_rule, \
            "sample_name lacks the genome suffix ShortStack's columns carry"

    def test_missing_columns_raise_a_named_error(self, srna_prep_rule):
        assert "missing" in srna_prep_rule and "raise ValueError" in srna_prep_rule, (
            "a column mismatch must name the missing columns, not surface as "
            "a bare pandas KeyError")

    def test_renames_strip_shortstack_decoration(self, srna_prep_rule):
        """The renames must leave exactly the mapped_name behind."""
        col = "clean__WT_sRNA_rep1__Spombe_condensed"
        assert col.removeprefix("clean__").removesuffix("_condensed") == \
            "WT_sRNA_rep1__Spombe"
        assert 'x[7:] if x.startswith("clean__")' in srna_prep_rule
        assert 'x[:-10] if x.endswith("_condensed")' in srna_prep_rule
