"""A handful of '.' rows must not cost a whole annotation its strandedness.

The `is_stranded` checkpoint decides whether a target file gets split by strand.
That split is what lets `define_key_for_plots` pick a sense-strand bigwig per
half, so demoting a file to "unstranded" silently turns stranded RNA and sRNA
tracks into separate plus/minus columns.

It used to demote on any value outside {'+', '-'}: 631 rows out of ColCEN's
42,927 TEs -- 1.5%, mostly GFF3's '?' for "stranded but unknown" -- were enough
to do it, while the gene annotation next to it was unaffected because it happens
to be pure ±.

The checkpoint body is exec'd out of the .smk rather than reimplemented here.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"


def _classify(tmp_path, strands, header=False):
    """Run the checkpoint's own run: body over a BED with these strand values."""
    text = SMK.read_text()
    start = text.index("checkpoint is_stranded:")
    block = text[start:text.index("###", start)]
    body = block[block.index("    run:") + len("    run:"):]
    body = "\n".join(l[8:] if l.startswith(" " * 8) else l
                     for l in body.splitlines())

    bed = tmp_path / "t.bed"
    rows = [f"chr1\t{i*10}\t{i*10+5}\tr{i}\t0\t{s}" for i, s in enumerate(strands)]
    if header:
        rows.insert(0, "#chrom\tstart\tend\tname\tscore\tstrand")
    bed.write_text("\n".join(rows) + "\n" if rows else "")
    hdr = tmp_path / "t.bed.header"
    hdr.write_text("yes\n" if header else "no\n")
    out = tmp_path / "t.bed.stranded"

    ns = {
        "input": type("I", (), {"bedfile": str(bed), "header": str(hdr)})(),
        "output": type("O", (), {"file": str(out)})(),
    }
    exec(compile(body, str(SMK), "exec"), ns)
    return out.read_text().strip()


class TestIsStranded:
    def test_pure_plus_minus_is_stranded(self, tmp_path):
        assert _classify(tmp_path, ["+", "-", "+", "-"]) == "stranded"

    def test_a_few_dots_stay_stranded(self, tmp_path):
        """The regression: 1.5% unstranded rows used to demote the whole file."""
        strands = ["+"] * 60 + ["-"] * 30 + ["."] * 4 + ["?"] * 2
        assert _classify(tmp_path, strands) == "stranded_mixed"

    def test_question_marks_alone_still_count_as_mixed(self, tmp_path):
        assert _classify(tmp_path, ["+", "-", "?"]) == "stranded_mixed"

    def test_no_oriented_rows_is_unstranded(self, tmp_path):
        assert _classify(tmp_path, [".", ".", "?"]) == "unstranded"

    def test_all_dots_is_unstranded(self, tmp_path):
        assert _classify(tmp_path, ["."] * 10) == "unstranded"

    def test_empty_file_is_unstranded(self, tmp_path):
        assert _classify(tmp_path, []) == "unstranded"

    @pytest.mark.parametrize("header", [False, True])
    def test_verdict_is_one_of_three(self, header, tmp_path):
        verdict = _classify(tmp_path, ["+", "-", "."], header=header)
        assert verdict in {"stranded", "stranded_mixed", "unstranded"}

    def test_header_row_is_not_read_as_a_strand(self, tmp_path):
        """With the header skipped, a pure ± file stays plain 'stranded'."""
        assert _classify(tmp_path, ["+", "-", "+"], header=True) == "stranded"
