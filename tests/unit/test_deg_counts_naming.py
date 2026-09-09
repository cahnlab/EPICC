"""prep_files_for_DEGs must read per-replicate counts by mapped_name.

The counts files are a post-alignment product, so they carry the genome. The
rule is a `run:` block, so a dry-run never executes it and the mismatch only
surfaced during a real run (#71).
"""

import os
import re
import sys

import pytest

_REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_RNASEQ_SMK = os.path.join(_REPO_ROOT, "workflow", "rules", "RNAseq.smk")

sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))
from scripts.sample_sheet import (  # noqa: E402
    add_compat_columns, read_sample_sheet,
)


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
