"""Tests for workflow/scripts/column_for_label.awk.

The header format is deeptools' multiBigwigSummary --outRawCounts output:
every name single-quoted, the first prefixed with '#'.
"""

import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "scripts" / "column_for_label.awk"

# One label is a substring of two others, which a regex match cannot tell apart.
HEADER = "#'chr'\t'start'\t'end'\t'WT_H3K9me2'\t'WT_H3K9me'\t'mut_WT_H3K9me'\t'WT_RNAseq_plus'\t'WT_RNAseq_minus'\n"
ROW = "chr1\t0\t100\t1.0\t2.0\t3.0\t4.0\t5.0\n"


def _lookup(tmp_path, label):
    tab = tmp_path / "values.tab"
    tab.write_text(HEADER + ROW)
    return subprocess.run(["awk", "-v", f"label={label}", "-f", str(SCRIPT), str(tab)],
                          capture_output=True, text=True)


@pytest.mark.parametrize("label,col", [
    ("WT_H3K9me2", "4"),
    ("WT_H3K9me", "5"),
    ("mut_WT_H3K9me", "6"),
    ("WT_RNAseq_plus", "7"),
    ("WT_RNAseq_minus", "8"),
])
def test_exact_match(tmp_path, label, col):
    result = _lookup(tmp_path, label)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == col


def test_missing_label_fails(tmp_path):
    result = _lookup(tmp_path, "WT_H3K27me3")
    assert result.returncode == 1
    assert "matches 0 columns" in result.stderr
    assert result.stdout == ""


def test_duplicate_label_fails(tmp_path):
    tab = tmp_path / "values.tab"
    tab.write_text("#'chr'\t'start'\t'end'\t'WT'\t'WT'\n" + "chr1\t0\t1\t1\t2\n")
    result = subprocess.run(["awk", "-v", "label=WT", "-f", str(SCRIPT), str(tab)],
                            capture_output=True, text=True)
    assert result.returncode == 1
    assert "matches 2 columns" in result.stderr


def test_prefix_and_quotes_stripped(tmp_path):
    # The first name carries both the '#' and the quotes.
    result = _lookup(tmp_path, "chr")
    assert result.stdout.strip() == "1"
