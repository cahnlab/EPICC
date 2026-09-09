"""Genome as a separate token: multi-genome sheets and mapped_name (issue #39).

Sample_ID identifies the library, not the alignment. Read processing stays
genome-free so one download+trim serves every reference; everything from
alignment onward carries the genome in `mapped_name`.
"""

import os
import sys

import pandas as pd
import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))

from sample_sheet import (  # noqa: E402
    add_compat_columns,
    explode_genomes,
    get_replicate_sample_ids,
    identify_control_samples,
    parse_genomes,
    read_sample_sheet,
)
from scripts.samplefile_validation import check_table  # noqa: E402

HEADER = ("Sample_ID\tAssay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
          "Read_layout\tIP_target\tControl")


def _sheet(tmp_path, rows, name="s.tsv"):
    f = tmp_path / name
    f.write_text("\n".join([HEADER] + rows) + "\n")
    return f


class TestParseGenomes:
    @pytest.mark.parametrize("raw,expected", [
        ("ColCEN", ["ColCEN"]),
        ("B73,W22", ["B73", "W22"]),
        (" B73 , W22 ", ["B73", "W22"]),
        ("B73,,W22", ["B73", "W22"]),
    ])
    def test_splits(self, raw, expected):
        assert parse_genomes(raw) == expected


class TestExplode:
    def test_one_row_per_genome(self):
        df = pd.DataFrame([{"Sample_ID": "s1", "Genome": "B73,W22"}])
        out = explode_genomes(df)
        assert list(out["Genome"]) == ["B73", "W22"]
        assert list(out["Sample_ID"]) == ["s1", "s1"]

    def test_single_genome_untouched(self):
        df = pd.DataFrame([{"Sample_ID": "s1", "Genome": "B73"}])
        assert list(explode_genomes(df)["Genome"]) == ["B73"]

    def test_repeated_reference_rejected(self):
        # After the explode a repeat is indistinguishable from two real rows,
        # so it has to be caught here.
        df = pd.DataFrame([{"Sample_ID": "s1", "Genome": "B73,B73"}])
        with pytest.raises(ValueError, match="more than once"):
            explode_genomes(df)


class TestMappedName:
    def test_mapped_name_is_sample_and_genome(self, tmp_path):
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73,W22\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        df = add_compat_columns(read_sample_sheet(f))
        assert sorted(df["mapped_name"]) == ["s1__B73", "s1__W22"]
        # sample_name stays bare: it is the pre-alignment identity.
        assert set(df["sample_name"]) == {"s1"}

    def test_single_genome_mapped_name(self, tmp_path):
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        df = add_compat_columns(read_sample_sheet(f))
        assert list(df["mapped_name"]) == ["s1__B73"]


class TestValidation:
    def test_multi_genome_sheet_validates(self, tmp_path):
        # Same Sample_ID and same reads on both exploded rows is the whole point.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73,W22\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        check_table(read_sample_sheet(f), check_paths=False)  # no raise

    def test_genome_with_double_underscore_rejected(self, tmp_path):
        # '__' is the delimiter in '{Sample_ID}__{Genome}'.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73__v5\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        with pytest.raises(ValueError, match="must not contain '__'"):
            check_table(read_sample_sheet(f), check_paths=False)

    def test_genuine_duplicate_still_rejected(self, tmp_path):
        # Same Sample_ID AND same genome is still a real duplicate.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
            "s1\tRNAseq\tB73\tgenotype:WT\trep2\tSRR2\tSE\t\t",
        ])
        with pytest.raises(ValueError, match="duplicate Sample_ID"):
            check_table(read_sample_sheet(f), check_paths=False)

    def test_shared_reads_across_samples_still_rejected(self, tmp_path):
        # Two DIFFERENT samples pointing at one accession remains an error.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
            "s2\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        with pytest.raises(ValueError, match="is also used by"):
            check_table(read_sample_sheet(f), check_paths=False)


class TestControlMergeAcceptsQualifiedNames:
    """A control merge group is keyed on a mapped_name, not a bare Sample_ID.

    An empty input list is legal in Snakemake, so a failed lookup only shows up
    as `samtools merge` with no inputs at runtime (#71) — dry-runs pass.
    """

    ROWS = [
        "ip1\tChIP\tB73,W22\tgenotype:WT\trep1\tSRR0000001\tSE\tH3K9me2\tin1",
        "ip2\tChIP\tB73,W22\tgenotype:WT\trep2\tSRR0000002\tSE\tH3K9me2\tin2",
        "in1\tChIP\tB73,W22\tgenotype:WT\trep1\tSRR0000003\tSE\tInput\t",
        "in2\tChIP\tB73,W22\tgenotype:WT\trep2\tSRR0000004\tSE\tInput\t",
    ]

    def _df(self, tmp_path):
        return add_compat_columns(read_sample_sheet(_sheet(tmp_path, self.ROWS)))

    def test_qualified_control_name_resolves(self, tmp_path):
        df = self._df(tmp_path)
        assert get_replicate_sample_ids("in1__B73", df) == ["in1__B73", "in2__B73"]

    def test_bare_control_name_still_resolves(self, tmp_path):
        df = self._df(tmp_path)
        assert get_replicate_sample_ids("in1", df) != []

    def test_replicates_never_mix_genomes(self, tmp_path):
        df = self._df(tmp_path)
        for genome in ("B73", "W22"):
            got = get_replicate_sample_ids("in1__" + genome, df)
            assert got, f"no replicates for in1__{genome}"
            assert {g.rsplit("__", 1)[1] for g in got} == {genome}

    def test_either_replicate_names_the_same_group(self, tmp_path):
        df = self._df(tmp_path)
        assert (get_replicate_sample_ids("in1__B73", df)
                == get_replicate_sample_ids("in2__B73", df))

    def test_unknown_name_still_returns_empty(self, tmp_path):
        df = self._df(tmp_path)
        assert get_replicate_sample_ids("nope__B73", df) == []


class TestEveryControlInEveryFixtureResolves:
    """No control may resolve to an empty replicate list, in any fixture."""

    SHEETS = [
        "test_samples_pombe.tsv",
        "test_samples_colcen.tsv",
        "test_samples_hg38_chr21.tsv",
    ]

    @pytest.mark.parametrize("sheet", SHEETS)
    def test_controls_resolve_to_their_replicates(self, sheet):
        path = os.path.join(_REPO_ROOT, "tests", "integration", "data", sheet)
        if not os.path.exists(path):
            pytest.skip(f"{sheet} not present")
        df = add_compat_columns(read_sample_sheet(path))
        controls = set(identify_control_samples(df))
        if not controls:
            pytest.skip(f"{sheet} declares no controls")
        checked = 0
        for _, row in df[df["Sample_ID"].isin(controls)].iterrows():
            name = row["mapped_name"]
            got = get_replicate_sample_ids(name, df)
            assert got, f"{sheet}: control {name} resolved to no replicates"
            assert name in got, f"{sheet}: {name} missing from its own group {got}"
            checked += 1
        assert checked > 0
