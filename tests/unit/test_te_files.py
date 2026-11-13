"""Tests for workflow/scripts/te_files.py: tagging an alternative TE annotation.

An alt_te_file (epicc --te-file) must never share a prepared file or an output
name with the options-file te_file, nor with a different alternative file.
"""

import hashlib
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "workflow"))

from scripts.te_files import (  # noqa: E402
    resolve_te_tag, te_bed_path, te_source, te_source_record,
    te_suffix_from_target, te_suffixes, te_tag_from_path, te_target_name,
)


@pytest.mark.parametrize("path,tag", [
    ("/data/EDTA_v2.gff3.gz", "EDTA_v2"),
    ("te.bed", "te"),
    ("https://host/x/Col-CEN.TEanno.gff3.gz?raw=true", "Col-CEN"),
    ("/d/my TEs (v3).bed", "my-TEs-v3"),
    ("/d/a__b.bed", "a_b"),          # '__' is the field separator
    ("/d/___.bed", "alt"),
])
def test_tag_from_path(path, tag):
    assert te_tag_from_path(path) == tag


def test_suffix_empty_without_alt(tmp_path):
    cfg = {"G": {"te_file": "a.bed"}, "H": {}}
    assert te_suffixes(cfg, str(tmp_path)) == {"G": "", "H": ""}


def test_suffix_from_alt(tmp_path):
    cfg = {"G": {"te_file": "a.bed", "alt_te_file": "/x/EDTA.bed"}}
    assert te_suffixes(cfg, str(tmp_path)) == {"G": "_EDTA"}


def test_alt_equal_to_te_file_is_untagged(tmp_path):
    cfg = {"G": {"te_file": "/x/a.bed", "alt_te_file": "/x/a.bed"}}
    assert te_suffixes(cfg, str(tmp_path)) == {"G": ""}


def test_same_source_keeps_tag(tmp_path):
    rec = Path(te_source_record(str(tmp_path), "G", "EDTA"))
    rec.parent.mkdir(parents=True)
    rec.write_text("/x/EDTA.bed\n")
    assert resolve_te_tag("/x/EDTA.bed", str(tmp_path), "G") == "EDTA"


def test_different_source_same_name_gets_hash(tmp_path):
    rec = Path(te_source_record(str(tmp_path), "G", "EDTA"))
    rec.parent.mkdir(parents=True)
    rec.write_text("/old/EDTA.bed\n")
    h = hashlib.sha1(b"/new/EDTA.bed").hexdigest()[:8]
    assert resolve_te_tag("/new/EDTA.bed", str(tmp_path), "G") == f"EDTA-{h}"


def test_source_and_paths():
    g = {"te_file": "a.bed", "alt_te_file": "b.bed"}
    assert te_source(g, "") == "a.bed"
    assert te_source(g, "_b") == "b.bed"
    assert te_bed_path("genomes", "G", "") == "genomes/G/G__TE_file.bed"
    assert te_bed_path("genomes", "G", "_b") == "genomes/G/G__TE_file_b.bed"
    assert te_target_name("") == "all_TEs"
    assert te_target_name("_b") == "all_TEs_b"


@pytest.mark.parametrize("target,suffix", [
    ("all_TEs", ""),
    ("all_TEs_EDTA_v2", "_EDTA_v2"),
    ("all_TEs_EDTA-3fa1c2d9", "_EDTA-3fa1c2d9"),
    ("all_TEsX", None),
    ("all_TEs_", None),
    ("all_genes", None),
])
def test_suffix_from_target(target, suffix):
    assert te_suffix_from_target(target) == suffix


def test_target_round_trip(tmp_path):
    cfg = {"G": {"alt_te_file": "/x/EDTA_v2.gff3.gz"}}
    suffix = te_suffixes(cfg, str(tmp_path))["G"]
    assert te_suffix_from_target(te_target_name(suffix)) == suffix


def test_shared_te_shell_unchanged_for_default_rule():
    # check_te_file's shell must stay byte-identical, or Snakemake's code
    # rerun-trigger re-prepares every TE file and all TE analyses after it.
    smk = (_REPO_ROOT / "workflow/rules/environment_setup.smk").read_text()
    rule = smk[smk.index("rule check_te_file:"):smk.index("rule check_alt_te_file:")]
    assert "shell:\n        TE_PREP_SHELL\n" in rule
    alt = smk[smk.index("rule check_alt_te_file:"):]
    assert "TE_PREP_SHELL +" in alt.split("\nrule ")[0]
