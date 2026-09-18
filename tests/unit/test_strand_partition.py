"""The plus/minus/nostrand passes must partition a target file exactly.

`making_stranded_matrix_on_targetfile` runs computeMatrix once per strand and
`merging_matrix` rbinds the results, so every region has to land in exactly one
pass: a row counted twice is plotted twice, and a row counted zero times
disappears from the figure without a word.

Annotations routinely carry a few rows that are neither '+' nor '-' -- GFF3 uses
'?' for "stranded but unknown", and 631 of ColCEN's 42,927 TEs are '.' or '?'.
Those go to the 'nostrand' pass. Note its predicate is the only one of the three
that a header line also satisfies, which is why it carries an explicit guard.

The awk text is taken out of the .smk rather than copied, so an edit there is
what these tests actually check.
"""

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"

STRANDS = ("plus", "minus", "nostrand", "unstranded")

BED_ROWS = [
    ("chr1", "100", "200", "plus_a", "0", "+"),
    ("chr1", "300", "400", "minus_a", "0", "-"),
    ("chr1", "500", "600", "dot_a", "0", "."),
    ("chr1", "700", "800", "question_a", "0", "?"),
    ("chr1", "900", "1000", "plus_b", "0", "+"),
    ("chr1", "1100", "1200", "minus_b", "0", "-"),
    ("chr1", "1300", "1400", "dot_b", "0", "."),
]
HEADER_ROW = ("#chrom", "start", "end", "name", "score", "strand")


def _render_filter(strand, target_file, out_file):
    """Extract the strand-splitting head of the rule and expand it like Snakemake.

    Everything up to the computeMatrix call: reading the header flag, the
    case/awk partition, and the '?' -> '.' normalisation. The computeMatrix call
    itself needs bigwigs, so it is left out.
    """
    text = SMK.read_text()
    start = text.index("rule making_stranded_matrix_on_targetfile:")
    block = text[start:text.index("rule merging_matrix:", start)]
    body = re.search(r'shell:\s*"""(.*?)"""', block, re.DOTALL).group(1)
    body = body[: body.index('echo "{params.labels}"')]
    body = body[body.index('header="$(cat'):]

    for placeholder, value in (
        ("{input.header}", str(target_file) + ".header"),
        ("{input.target_file}", str(target_file)),
        ("{output.temp}", str(out_file)),
        ("{params.strand}", strand),
    ):
        body = body.replace(placeholder, value)
    return body.replace("{{", "{").replace("}}", "}")


def _run(strand, tmp_path, rows, header=False):
    bed = tmp_path / "target.bed"
    lines = ([HEADER_ROW] if header else []) + list(rows)
    bed.write_text("\n".join("\t".join(r) for r in lines) + "\n")
    (tmp_path / "target.bed.header").write_text("yes\n" if header else "no\n")

    out = tmp_path / f"temp_{strand}.bed"
    script = _render_filter(strand, bed, out)
    workdir = tmp_path / "tmpdir"
    workdir.mkdir(exist_ok=True)
    proc = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True,
        env={"PATH": __import__("os").environ["PATH"], "TMPDIR": str(workdir)},
    )
    assert proc.returncode == 0, f"{strand} pass failed:\n{proc.stdout}\n{proc.stderr}"
    text = out.read_text()
    return [l.split("\t") for l in text.splitlines() if l]


def _names(rows):
    return [r[3] for r in rows]


class TestStrandPartition:
    def test_plus_pass_takes_only_plus(self, tmp_path):
        assert _names(_run("plus", tmp_path, BED_ROWS)) == ["plus_a", "plus_b"]

    def test_minus_pass_takes_only_minus(self, tmp_path):
        assert _names(_run("minus", tmp_path, BED_ROWS)) == ["minus_a", "minus_b"]

    def test_nostrand_pass_takes_dot_and_question(self, tmp_path):
        assert _names(_run("nostrand", tmp_path, BED_ROWS)) == [
            "dot_a", "question_a", "dot_b"
        ]

    def test_three_passes_partition_the_file(self, tmp_path):
        """Every region in exactly one pass -- none lost, none duplicated."""
        seen = []
        for strand in ("plus", "minus", "nostrand"):
            seen += _names(_run(strand, tmp_path, BED_ROWS))
        assert sorted(seen) == sorted(r[3] for r in BED_ROWS)
        assert len(seen) == len(set(seen))

    def test_unstranded_pass_takes_everything(self, tmp_path):
        assert _names(_run("unstranded", tmp_path, BED_ROWS)) == [
            r[3] for r in BED_ROWS
        ]

    @pytest.mark.parametrize("strand", STRANDS)
    def test_header_row_is_never_emitted(self, strand, tmp_path):
        """'nostrand' is the one predicate a header line also matches."""
        rows = _run(strand, tmp_path, BED_ROWS, header=True)
        assert HEADER_ROW[0] not in [r[0] for r in rows]

    @pytest.mark.parametrize("strand", STRANDS)
    def test_strand_column_is_normalised(self, strand, tmp_path):
        """computeMatrix panics on '?', so nothing may reach it unnormalised."""
        for row in _run(strand, tmp_path, BED_ROWS):
            assert row[5] in {"+", "-", "."}, f"{strand}: bad strand {row[5]!r}"

    def test_all_unstranded_file_yields_nothing_for_plus_and_minus(self, tmp_path):
        """A file with no ± at all is classified unstranded and never split."""
        rows = [r for r in BED_ROWS if r[5] not in {"+", "-"}]
        assert _run("plus", tmp_path, rows) == []
        assert _run("minus", tmp_path, rows) == []
        assert _names(_run("nostrand", tmp_path, rows)) == [r[3] for r in rows]
