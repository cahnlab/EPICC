"""Tests for the methylation_contexts check in bedmethyl_to_cx_report.py.

A dmC sample whose bedMethyl has no covered call in a context listed in
methylation_contexts (a CpG-only ONT model under CG/CHG/CHH) must stop at
the conversion, naming methylation_contexts, rather than at DMR calling.
"""

import sys
import types
from pathlib import Path

import pytest


class _FastaFile:
    """Minimal stand-in for pysam.FastaFile (fetch + close) on a plain FASTA."""

    def __init__(self, path):
        self.seqs, name = {}, None
        for line in Path(path).read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                self.seqs[name] = ""
            elif name:
                self.seqs[name] += line.strip()

    def fetch(self, chrom):
        return self.seqs[chrom]

    def close(self):
        pass


# The script only needs FastaFile; the epicc env that runs unit tests has no
# pysam, so stand one in when it is missing.
try:
    import pysam  # noqa: F401
    _stubbed = False
except ImportError:
    sys.modules["pysam"] = types.SimpleNamespace(FastaFile=_FastaFile)
    _stubbed = True

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))

from bedmethyl_to_cx_report import (  # noqa: E402
    convert_bedmethyl_to_cx_report,
    missing_contexts_error,
)

# The module keeps its own reference; don't let the stub reach other tests.
if _stubbed:
    del sys.modules["pysam"]

# 0-based: C at 2 is CG, C at 6 is CHG (CAG), C at 10 is CHH (CAT).
SEQ = "AACGTTCAGTTCATTAA"
CG_POS, CHG_POS, CHH_POS = 2, 6, 10


def _row(pos, n_mod, n_canon, code="m"):
    cov = n_mod + n_canon
    # chrom start end code score strand tstart tend color Nvalid frac Nmod Ncanon ...
    return "\t".join(map(str, [
        "chr1", pos, pos + 1, code, cov, "+", pos, pos + 1, "255,0,0",
        cov, round(100 * n_mod / cov, 2) if cov else 0, n_mod, n_canon, 0, 0, 0, 0, 0,
    ]))


@pytest.fixture
def fasta(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{SEQ}\n")
    return str(fa)


def _bedmethyl(tmp_path, rows):
    p = tmp_path / "in.bedmethyl"
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def test_cpg_only_with_all_contexts_required_stops(tmp_path, fasta, capsys):
    bm = _bedmethyl(tmp_path, [_row(CG_POS, 5, 5)])
    with pytest.raises(SystemExit) as exc:
        convert_bedmethyl_to_cx_report(bm, fasta, str(tmp_path / "out.cx"),
                                       require_contexts=["CG", "CHG", "CHH"],
                                       sample="s1")
    assert exc.value.code == 1
    err = capsys.readouterr().err
    assert "dmC sample s1 has no covered CHG/CHH position" in err
    assert "methylation_contexts" in err
    assert "all-context model" in err


def test_cpg_only_with_cg_required_passes(tmp_path, fasta):
    bm = _bedmethyl(tmp_path, [_row(CG_POS, 5, 5)])
    out = tmp_path / "out.cx"
    covered = convert_bedmethyl_to_cx_report(bm, fasta, str(out),
                                             require_contexts=["CG"])
    assert covered == {"CG": 1, "CHG": 0, "CHH": 0}
    assert out.read_text().split("\t")[5] == "CG"


def test_no_requirement_never_stops(tmp_path, fasta):
    bm = _bedmethyl(tmp_path, [_row(CG_POS, 5, 5)])
    covered = convert_bedmethyl_to_cx_report(bm, fasta, str(tmp_path / "o.cx"))
    assert covered["CG"] == 1


def test_all_context_data_passes(tmp_path, fasta):
    bm = _bedmethyl(tmp_path, [_row(CG_POS, 5, 5), _row(CHG_POS, 1, 3),
                               _row(CHH_POS, 0, 4)])
    covered = convert_bedmethyl_to_cx_report(
        bm, fasta, str(tmp_path / "o.cx"), require_contexts=["CG", "CHG", "CHH"])
    assert covered == {"CG": 1, "CHG": 1, "CHH": 1}


def test_zero_coverage_row_does_not_count(tmp_path, fasta, capsys):
    bm = _bedmethyl(tmp_path, [_row(CG_POS, 5, 5), _row(CHG_POS, 0, 0)])
    with pytest.raises(SystemExit):
        convert_bedmethyl_to_cx_report(bm, fasta, str(tmp_path / "o.cx"),
                                       require_contexts=["CG", "CHG"], sample="s2")
    assert "no covered CHG position, but methylation_contexts lists it" in capsys.readouterr().err


def test_5hmc_rows_do_not_count(tmp_path, fasta):
    bm = _bedmethyl(tmp_path, [_row(CG_POS, 5, 5), _row(CHG_POS, 3, 3, code="h")])
    with pytest.raises(SystemExit):
        convert_bedmethyl_to_cx_report(bm, fasta, str(tmp_path / "o.cx"),
                                       require_contexts=["CG", "CHG"])


def test_error_lists_contexts_with_data():
    msg = missing_contexts_error("s", ["CHH"], {"CG": 10, "CHG": 2, "CHH": 0})
    assert "(calls per context: CG 10, CHG 2)" in msg
