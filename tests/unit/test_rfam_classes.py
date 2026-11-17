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


def depleted_names(filt, keep=()):
    bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(), filt, keep)
    return [line.split("\t")[3].split("::")[0] for line in bed]


def cfg(**opt):
    return {"structural_rna_classes": opt}


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


class TestDefaultSelection:
    def test_default_is_rrna_trna_snrna(self):
        assert set(rc.DEFAULT_FILTER) == {"rRNA", "tRNA", "snRNA"}
        assert rc.DEFAULT_KEEP == ()
        assert depleted_names(rc.DEFAULT_FILTER) == [
            "tRNA", "snoR71", "U1", "5S_rRNA"]

    def test_default_takes_all_snrnas_and_keeps_mirna(self):
        """snRNA is snoRNA's parent: the default takes spliceosomal snRNAs and
        every snoRNA subtype, and leaves miRNA alone."""
        fam = families()
        for acc in ("RF00003",            # Gene; snRNA; splicing;
                    "RF00097", "RF00012",  # CD-box
                    "RF00090",            # HACA-box
                    "RF00286"):           # scaRNA
            assert rc.is_depleted(fam[acc][1], rc.DEFAULT_FILTER), acc
        assert not rc.is_depleted(fam["RF00103"][1], rc.DEFAULT_FILTER)

    def test_default_spares_mirna_cisreg_and_unknown(self):
        names = depleted_names(rc.DEFAULT_FILTER)
        for spared in ("mir-1", "Cobalamin", "new_fam"):
            assert spared not in names

    def test_snorna_matches_all_snorna_subtypes_only(self):
        fam = families()
        for acc in ("RF00097", "RF00090", "RF00286"):
            assert rc.is_depleted(fam[acc][1], {"snoRNA"})
        assert not rc.is_depleted(fam["RF00003"][1], {"snoRNA"})

    def test_parent_term_selects_whole_branch(self):
        assert depleted_names({"snRNA"}) == ["snoR71", "U1"]


class TestFilterAndKeep:
    """A family is depleted when a term is in filter and none is in keep."""

    def test_keep_exempts_a_subtree(self):
        # snRNA minus the snoRNAs leaves only the spliceosomal U1
        assert depleted_names({"snRNA"}, {"snoRNA"}) == ["U1"]

    def test_family_with_terms_in_both_lists_is_kept(self):
        fam = families()
        terms = fam["RF00097"][1]  # Gene; snRNA; snoRNA; CD-box
        assert rc.is_depleted(terms, {"snRNA"})
        assert not rc.is_depleted(terms, {"snRNA"}, {"CD-box"})
        assert not rc.is_depleted(terms, {"CD-box"}, {"snRNA"})  # keep wins either way

    def test_everything_but_mirna(self):
        names = depleted_names({"Gene", "Cis-reg", "Intron"}, {"miRNA"})
        assert "mir-1" not in names
        assert set(names) == {"5S_rRNA", "tRNA", "snoR71", "U1", "Cobalamin"}

    def test_unclassified_accession_never_depleted(self):
        assert "new_fam" not in depleted_names({"Gene", "Cis-reg", "Intron"})
        assert not rc.is_depleted(None, {"Gene"})

    def test_empty_filter_depletes_nothing(self):
        assert depleted_names(()) == []
        assert depleted_names((), {"miRNA"}) == []

    def test_keep_alone_applies_to_default_filter(self):
        assert rc.selection_from_config(cfg(keep=["snoRNA"])) == (
            ("rRNA", "snRNA", "tRNA"), ("snoRNA",))

    def test_tallies_and_missing(self):
        _, dep, spared, missing = rc.select_hits(
            TBLOUT.splitlines(), families(), rc.DEFAULT_FILTER)
        assert dep["Gene; rRNA"] == 1 and dep["Gene; tRNA"] == 1
        assert spared["Gene; miRNA"] == 1
        assert spared["(not in family table)"] == 1
        assert missing == {"RF09999"}


class TestOverlap:
    def test_keep_wins_and_warns(self):
        c = cfg(filter=["rRNA", "tRNA", "snRNA"], keep=["tRNA"])
        assert rc.selection_from_config(c) == (("rRNA", "snRNA"), ("tRNA",))
        warns = rc.structural_rna_warnings(c)
        assert len(warns) == 1
        assert "tRNA" in warns[0] and "'keep' wins" in warns[0]
        assert rc.validate_structural_rna_classes(c) == []

    def test_overlap_of_whole_filter_disables_depletion(self):
        c = cfg(filter=["snRNA"], keep=["snRNA"])
        assert rc.selection_from_config(c) == ((), ("snRNA",))
        assert not rc.depletion_enabled(c)

    def test_no_warning_without_overlap(self):
        assert rc.structural_rna_warnings(cfg(filter=["rRNA"], keep=["miRNA"])) == []
        assert rc.structural_rna_warnings({}) == []


class TestEnabled:
    def test_default_enabled(self):
        assert rc.depletion_enabled({})
        assert rc.depletion_enabled(cfg(keep=["miRNA"]))

    def test_empty_filter_disables(self):
        assert not rc.depletion_enabled(cfg(filter=[]))
        assert not rc.depletion_enabled(cfg(filter=[], keep=["miRNA"]))


class TestDeprecatedToggle:
    def test_false_means_empty_filter(self):
        c = {"structural_rna_depletion": False}
        assert not rc.depletion_enabled(c)
        assert rc.selection_from_config(c)[0] == ()
        warns = rc.structural_rna_warnings(c)
        assert len(warns) == 1 and "deprecated" in warns[0]
        assert "filter: []" in warns[0]

    def test_false_overrides_shipped_filter(self):
        # the shipped options file always merges in a filter list
        c = {"structural_rna_depletion": False,
             "structural_rna_classes": {"filter": ["rRNA", "tRNA", "snRNA"],
                                        "keep": []}}
        assert not rc.depletion_enabled(c)

    def test_true_is_ignored_with_warning(self):
        c = {"structural_rna_depletion": True}
        assert rc.selection_from_config(c) == rc.selection_from_config({})
        warns = rc.structural_rna_warnings(c)
        assert len(warns) == 1 and "ignored" in warns[0]
        assert not rc.depletion_enabled(
            {"structural_rna_depletion": True, **cfg(filter=[])})

    def test_non_bool_refused(self):
        errs = rc.validate_structural_rna_classes({"structural_rna_depletion": "no"})
        assert errs and "true or false" in errs[0]


class TestBedRecord:
    def test_minus_strand_coordinates_and_name(self):
        bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(), {"tRNA"})
        assert bed == ["chr1\t508\t580\ttRNA::chr1:508-580(-)\t0\t-"]

    def test_plus_strand(self):
        bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(), {"rRNA"})
        assert bed == ["chr2\t100\t220\t5S_rRNA::chr2:100-220(+)\t0\t+"]

    def test_sorted_by_chrom_then_start(self):
        bed, *_ = rc.select_hits(TBLOUT.splitlines(), families(),
                                 {"Gene", "Cis-reg"})
        keys = [(l.split("\t")[0], int(l.split("\t")[1])) for l in bed]
        assert keys == sorted(keys)


class TestSelectionLabel:
    def test_round_trip(self):
        for sel in [(("rRNA", "snRNA", "tRNA"), ()), ((), ("miRNA",)),
                    (("Gene",), ("miRNA", "snoRNA")), ((), ())]:
            label = rc.selection_label(*sel)
            assert rc.parse_selection_label(label) == sel

    def test_default_from_config(self):
        assert rc.selection_from_config({}) == (("rRNA", "snRNA", "tRNA"), ())

    def test_sorted_and_deduped(self):
        c = cfg(filter=["tRNA", "rRNA", "tRNA"], keep=["miRNA", "miRNA"])
        assert rc.selection_from_config(c) == (("rRNA", "tRNA"), ("miRNA",))

    def test_label_tracks_effective_selection_only(self):
        def label(c):
            return rc.selection_label(*rc.selection_from_config(c))
        assert label(cfg(filter=["tRNA", "rRNA"])) == label(cfg(filter=["rRNA", "tRNA"]))
        # omitted keep == empty keep; omitted filter == default filter
        assert label(cfg(filter=["rRNA"])) == label(cfg(filter=["rRNA"], keep=[]))
        assert label({}) == label(cfg(keep=[]))
        # a filter class that keep cancels contributes nothing
        assert label(cfg(filter=["rRNA", "tRNA"], keep=["tRNA"])) == \
            label(cfg(filter=["rRNA"], keep=["tRNA"]))
        assert label(cfg(filter=["rRNA"])) != label(cfg(filter=["rRNA"], keep=["miRNA"]))

    def test_bad_label(self):
        with pytest.raises(ValueError):
            rc.parse_selection_label("drop:rRNA")


class TestValidation:
    def v(self, opt):
        return rc.validate_structural_rna_classes({"structural_rna_classes": opt})

    def test_absent_or_null_is_valid(self):
        assert rc.validate_structural_rna_classes({}) == []
        assert self.v(None) == []
        assert self.v({}) == []

    def test_valid_combinations(self):
        assert self.v({"filter": ["rRNA", "tRNA", "snoRNA"]}) == []
        assert self.v({"keep": ["miRNA", "lncRNA"]}) == []
        assert self.v({"filter": [], "keep": []}) == []
        assert self.v({"filter": ["Gene"], "keep": ["miRNA"]}) == []
        assert self.v({"filter": ["rRNA"], "keep": ["rRNA"]}) == []  # warning only

    def test_unknown_class_in_filter(self):
        errs = self.v({"filter": ["rRNA", "snorna", "microRNA"]})
        assert len(errs) == 1
        assert "filter" in errs[0]
        assert "microRNA" in errs[0] and "snorna" in errs[0]
        assert "snoRNA" in errs[0]  # lists the known classes

    def test_unknown_class_in_keep(self):
        errs = self.v({"keep": ["miRNAs"]})
        assert len(errs) == 1 and "keep" in errs[0] and "miRNAs" in errs[0]

    def test_unknown_in_both_reported_separately(self):
        assert len(self.v({"filter": ["x"], "keep": ["y"]})) == 2

    def test_unknown_key(self):
        errs = self.v({"remove": ["rRNA"]})
        assert any("remove" in e for e in errs)

    def test_not_a_mapping(self):
        assert self.v(["rRNA"]) != []

    def test_string_or_null_instead_of_list(self):
        assert any("must be a list" in e for e in self.v({"filter": "rRNA"}))
        errs = self.v({"filter": None})
        assert any("[]" in e for e in errs)

    def test_known_classes_cover_defaults(self):
        assert set(rc.DEFAULT_FILTER) <= rc.KNOWN_CLASSES


def _shipped_options():
    import yaml
    path = os.path.join(os.path.dirname(__file__), "..", "..", "config",
                        "epicc-options.yaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


def test_shipped_options_match_code_defaults():
    shipped = _shipped_options()
    assert "structural_rna_depletion" not in shipped
    assert rc.validate_structural_rna_classes(shipped) == []
    assert rc.structural_rna_warnings(shipped) == []
    assert rc.selection_from_config(shipped) == rc.selection_from_config({})


def test_user_options_deep_merge_over_shipped():
    """Snakemake deep-merges a user --configfile over the shipped options:
    a user keep combines with the shipped filter, a user list replaces."""
    update_config = pytest.importorskip("snakemake.utils").update_config
    import copy

    def merged(user):
        c = copy.deepcopy(_shipped_options())
        update_config(c, user)
        return c

    assert rc.selection_from_config(merged(cfg(keep=["snoRNA"]))) == (
        ("rRNA", "snRNA", "tRNA"), ("snoRNA",))
    assert not rc.depletion_enabled(merged(cfg(filter=[])))
    assert not rc.depletion_enabled(merged({"structural_rna_depletion": False}))


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
         "--selection", "filter:rRNA,snRNA,tRNA;keep:snoRNA",
         "--out-bed", str(bed)],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    names = [l.split("\t")[3].split("::")[0] for l in bed.read_text().splitlines()]
    assert names == ["tRNA", "U1", "5S_rRNA"]
    assert "Depleting 3 hits; sparing 4" in res.stderr
    assert "RF09999" in res.stderr  # warned about the unknown accession
