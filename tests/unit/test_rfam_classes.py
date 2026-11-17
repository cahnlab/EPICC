"""Rfam class selection for sRNA structural-RNA depletion.

Rfam holds miRNA and other non-structural families alongside rRNA/tRNA/snRNA,
so the auto-derived depletion FASTA must be built from a chosen subset of
classes (rfam_classes.py, rule select_structural_rna_classes).
"""

import gzip
import os
import subprocess
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)

import rfam_classes as rc  # noqa: E402


def _family_row(acc, fam_id, type_str):
    """A 35-column family.txt row with `type` in column 19, as Rfam ships it."""
    fields = [""] * 35
    fields[0], fields[1], fields[3] = acc, fam_id, f"{fam_id} description"
    fields[9] = "free text comment with a stray \x92 byte"
    fields[18] = type_str
    return "\t".join(fields) + "\n"


FAMILY_ROWS = [
    ("RF00001", "5S_rRNA", "Gene; rRNA;"),
    ("RF00005", "tRNA", "Gene; tRNA;"),
    ("RF00097", "snoR71", "Gene; snRNA; snoRNA; CD-box;"),
    ("RF00012", "U3", "Gene; snRNA; snoRNA; CD-box;"),
    ("RF00090", "snoR10", "Gene; snRNA; snoRNA; HACA-box;"),
    ("RF00286", "scaRNA", "Gene; snRNA; snoRNA; scaRNA;"),
    ("RF00003", "U1", "Gene; snRNA; splicing;"),
    ("RF00103", "mir-1", "Gene; miRNA;"),
    ("RF00174", "Cobalamin", "Cis-reg; riboswitch;"),
    ("RF01960", "SSU_rRNA_eukarya", "Gene; rRNA;"),
    ("RF00017", "Metazoa_SRP", "Gene;"),
    ("RF02924", "antisense_x", "Gene; antisense; "),  # trailing space, as in Rfam
]
FAMILY_TEXT = "".join(_family_row(*r) for r in FAMILY_ROWS)


def _tbl(idx, fam, acc, chrom, start, end, strand):
    """One cmscan --fmt 2 tblout data line (description has spaces)."""
    return (f"{idx} {fam} {acc} {chrom} - - cm 1 100 {start} {end} {strand} "
            f"no 1 0.50 0.0 80.0 1e-20 ! * - - - - - - {fam} some description\n")


TBLOUT = (
    "#idx target name accession query name ...\n"
    "#--- -----------\n"
    + _tbl(1, "5S_rRNA", "RF00001", "chr2", 101, 220, "+")
    + _tbl(2, "tRNA", "RF00005", "chr1", 580, 509, "-")
    + _tbl(3, "snoR71", "RF00097", "chr1", 1001, 1080, "+")
    + _tbl(4, "U1", "RF00003", "chr1", 2001, 2160, "+")
    + _tbl(5, "mir-1", "RF00103", "chr1", 3001, 3090, "+")
    + _tbl(6, "Cobalamin", "RF00174", "chr3", 11, 200, "+")
    + _tbl(7, "new_fam", "RF09999", "chr3", 501, 600, "+")  # not in table
)


def families():
    return rc.read_family_table(FAMILY_TEXT.splitlines(keepends=True))


def depleted_names(mode, classes):
    bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(), mode, classes)
    return [line.split("\t")[3].split("::")[0] for line in bed]


class TestParseType:
    def test_splits_and_strips(self):
        assert rc.parse_type("Gene; snRNA; snoRNA; CD-box;") == (
            "Gene", "snRNA", "snoRNA", "CD-box")

    def test_trailing_whitespace_and_empty(self):
        assert rc.parse_type("Gene; antisense; ") == ("Gene", "antisense")
        assert rc.parse_type("") == ()
        assert rc.parse_type(None) == ()


class TestFamilyTable:
    def test_reads_type_column(self):
        fam = families()
        assert fam["RF00001"] == ("5S_rRNA", ("Gene", "rRNA"))
        assert fam["RF00097"][1] == ("Gene", "snRNA", "snoRNA", "CD-box")
        assert len(fam) == len(FAMILY_ROWS)

    def test_skips_short_and_non_family_lines(self):
        assert rc.read_family_table(["junk\n", "RF1\tx\n"]) == {}


class TestSelection:
    def test_default_is_rrna_trna_snrna(self):
        assert set(rc.DEFAULT_CLASSES) == {"rRNA", "tRNA", "snRNA"}
        assert depleted_names("filter", rc.DEFAULT_CLASSES) == [
            "tRNA", "snoR71", "U1", "5S_rRNA"]

    def test_default_takes_all_snrnas_and_keeps_mirna(self):
        """snRNA is snoRNA's parent: the default takes spliceosomal snRNAs and
        every snoRNA subtype, and leaves miRNA alone."""
        fam = families()
        for acc in ("RF00003",            # Gene; snRNA; splicing;
                    "RF00097", "RF00012",  # CD-box
                    "RF00090",            # HACA-box
                    "RF00286"):           # scaRNA
            assert rc.is_depleted(fam[acc][1], "filter", rc.DEFAULT_CLASSES), acc
        assert not rc.is_depleted(fam["RF00103"][1], "filter", rc.DEFAULT_CLASSES)

    def test_default_spares_mirna_cisreg_and_unknown(self):
        names = depleted_names("filter", rc.DEFAULT_CLASSES)
        for spared in ("mir-1", "Cobalamin", "new_fam"):
            assert spared not in names

    def test_snorna_matches_all_snorna_subtypes(self):
        fam = families()
        for acc in ("RF00097", "RF00090", "RF00286"):
            assert rc.is_depleted(fam[acc][1], "filter", {"snoRNA"})
        assert not rc.is_depleted(fam["RF00003"][1], "filter", {"snoRNA"})

    def test_parent_term_selects_whole_branch(self):
        names = depleted_names("filter", {"snRNA"})
        assert names == ["snoR71", "U1"]

    def test_keep_depletes_everything_else(self):
        names = depleted_names("keep", {"miRNA"})
        assert "mir-1" not in names
        assert set(names) == {"5S_rRNA", "tRNA", "snoR71", "U1", "Cobalamin",
                              "new_fam"}

    def test_keep_empty_depletes_all_hits(self):
        assert len(depleted_names("keep", ())) == 7

    def test_tallies_and_missing(self):
        _, dep, spared, missing = rc.select_hits(
            TBLOUT.splitlines(), families(), "filter", rc.DEFAULT_CLASSES)
        assert dep["Gene; rRNA"] == 1 and dep["Gene; tRNA"] == 1
        assert spared["Gene; miRNA"] == 1
        assert spared["(not in family table)"] == 1
        assert missing == {"RF09999"}


class TestBedRecord:
    def test_minus_strand_coordinates_and_name(self):
        bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(), "filter",
                                 {"tRNA"})
        assert bed == ["chr1\t508\t580\ttRNA::chr1:508-580(-)\t0\t-"]

    def test_plus_strand(self):
        bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(), "filter",
                                 {"rRNA"})
        assert bed == ["chr2\t100\t220\t5S_rRNA::chr2:100-220(+)\t0\t+"]

    def test_sorted_by_chrom_then_start(self):
        bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(), "keep", ())
        keys = [(l.split("\t")[0], int(l.split("\t")[1])) for l in bed]
        assert keys == sorted(keys)


class TestSelectionLabel:
    def test_round_trip(self):
        for mode, classes in [("filter", ("rRNA", "snoRNA", "tRNA")),
                              ("keep", ()), ("keep", ("miRNA",))]:
            label = rc.selection_label(mode, classes)
            assert rc.parse_selection_label(label) == (mode, classes)

    def test_default_from_config(self):
        assert rc.selection_from_config({}) == ("filter", ("rRNA", "snRNA", "tRNA"))

    def test_from_config_sorted_and_deduped(self):
        cfg = {"structural_rna_classes": {"filter": ["tRNA", "rRNA", "tRNA"]}}
        assert rc.selection_from_config(cfg) == ("filter", ("rRNA", "tRNA"))
        cfg = {"structural_rna_classes": {"keep": []}}
        assert rc.selection_from_config(cfg) == ("keep", ())

    def test_order_does_not_change_label(self):
        a = rc.selection_label(*rc.selection_from_config(
            {"structural_rna_classes": {"filter": ["tRNA", "rRNA"]}}))
        b = rc.selection_label(*rc.selection_from_config(
            {"structural_rna_classes": {"filter": ["rRNA", "tRNA"]}}))
        assert a == b

    def test_bad_mode(self):
        with pytest.raises(ValueError):
            rc.parse_selection_label("drop:rRNA")


class TestValidation:
    def v(self, opt):
        return rc.validate_structural_rna_classes({"structural_rna_classes": opt})

    def test_absent_or_null_is_valid(self):
        assert rc.validate_structural_rna_classes({}) == []
        assert self.v(None) == []

    def test_valid_filter_and_keep(self):
        assert self.v({"filter": ["rRNA", "tRNA", "snoRNA"]}) == []
        assert self.v({"keep": ["miRNA", "lncRNA"]}) == []
        assert self.v({"keep": []}) == []

    def test_unknown_class(self):
        errs = self.v({"filter": ["rRNA", "snorna", "microRNA"]})
        assert len(errs) == 1
        assert "microRNA" in errs[0] and "snorna" in errs[0]
        assert "snoRNA" in errs[0]  # lists the known classes

    def test_both_modes_refused(self):
        errs = self.v({"filter": ["rRNA"], "keep": ["miRNA"]})
        assert any("not both" in e for e in errs)

    def test_unknown_key(self):
        errs = self.v({"remove": ["rRNA"]})
        assert any("remove" in e for e in errs)

    def test_neither_mode(self):
        assert self.v({}) != []

    def test_not_a_mapping(self):
        assert self.v(["rRNA"]) != []

    def test_string_instead_of_list(self):
        errs = self.v({"filter": "rRNA"})
        assert any("must be a list" in e for e in errs)

    def test_empty_filter_refused(self):
        errs = self.v({"filter": []})
        assert any("structural_rna_depletion" in e for e in errs)

    def test_known_classes_cover_defaults(self):
        assert set(rc.DEFAULT_CLASSES) <= rc.KNOWN_CLASSES


def test_shipped_options_leave_the_option_commented():
    """The default must come from code, not the options file.

    Snakemake deep-merges a user --configfile over the shipped one, so a
    shipped `filter:` plus a user's `keep:` would merge into both keys and
    fail validation.
    """
    import yaml
    path = os.path.join(os.path.dirname(__file__), "..", "..", "config",
                        "epicc-options.yaml")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    assert rc.CONFIG_KEY not in cfg
    assert rc.validate_structural_rna_classes(cfg) == []


def test_cli_writes_bed(tmp_path):
    fam = tmp_path / "family.txt.gz"
    with gzip.open(fam, "wb") as fh:
        fh.write(FAMILY_TEXT.encode("latin-1"))  # raw 0x92, not UTF-8
    tbl = tmp_path / "infernal.tblout"
    tbl.write_text(TBLOUT)
    bed = tmp_path / "out.bed"
    res = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "rfam_classes.py"),
         "--family-table", str(fam), "--tblout", str(tbl),
         "--selection", "filter:rRNA,snRNA,tRNA", "--out-bed", str(bed)],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    names = [l.split("\t")[3].split("::")[0] for l in bed.read_text().splitlines()]
    assert names == ["tRNA", "snoR71", "U1", "5S_rRNA"]
    assert "Depleting 4 hits; sparing 3" in res.stderr
    assert "RF09999" in res.stderr  # warned about the unknown accession
