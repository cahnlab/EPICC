"""Tests for workflow/scripts/bedgraph_to_bigwig.sh.

Input grouped by chromosome with ascending starts is ordered by concatenating
per-chromosome blocks; anything else falls back to sort. Both paths must give
the same track, padded with a zero-value base on every chromosome the data
does not cover.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "scripts" / "bedgraph_to_bigwig.sh"

if not (shutil.which("bedGraphToBigWig") and shutil.which("bigWigToBedGraph")):
    pytest.skip("bedGraphToBigWig/bigWigToBedGraph not found", allow_module_level=True)

# chr10 sorts before chr2 in C order, so a file in this order is grouped but
# not sorted. chr3 has no data and must be padded.
CHROM_SIZES = "chr2\t1000\nchr10\t1000\nchr1\t1000\nchr3\t1000\n"

EXPECTED = [
    "chr1\t5\t6\t50",
    "chr1\t9\t10\t100",
    "chr10\t0\t1\t25",
    "chr10\t700\t701\t75",
    "chr2\t3\t4\t10",
    "chr2\t40\t41\t20",
    "chr3\t0\t1\t0",
]


def _convert(tmp_path, lines):
    bg = tmp_path / "in.bedGraph"
    bg.write_text("".join(f"{line}\n" for line in lines))
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text(CHROM_SIZES)
    bw = tmp_path / "out.bw"
    result = subprocess.run(
        ["bash", str(SCRIPT), str(bg), str(sizes), str(bw)],
        capture_output=True, text=True, env={"TMPDIR": str(tmp_path),
                                             "PATH": subprocess.os.environ["PATH"]},
    )
    assert result.returncode == 0, result.stderr
    back = subprocess.run(["bigWigToBedGraph", str(bw), "/dev/stdout"],
                          capture_output=True, text=True, check=True)
    return [line for line in back.stdout.splitlines() if line]


def test_grouped_input_in_genome_order(tmp_path):
    lines = [
        "chr2\t3\t4\t10", "chr2\t40\t41\t20",
        "chr10\t0\t1\t25", "chr10\t700\t701\t75",
        "chr1\t5\t6\t50", "chr1\t9\t10\t100",
    ]
    assert _convert(tmp_path, lines) == EXPECTED


def test_chromosome_revisited_falls_back_to_sort(tmp_path):
    lines = [
        "chr2\t3\t4\t10", "chr10\t0\t1\t25", "chr1\t5\t6\t50",
        "chr2\t40\t41\t20", "chr10\t700\t701\t75", "chr1\t9\t10\t100",
    ]
    assert _convert(tmp_path, lines) == EXPECTED


def test_descending_start_falls_back_to_sort(tmp_path):
    # 40 then 3 compares lower numerically but higher as a string.
    lines = [
        "chr2\t40\t41\t20", "chr2\t3\t4\t10",
        "chr10\t700\t701\t75", "chr10\t0\t1\t25",
        "chr1\t9\t10\t100", "chr1\t5\t6\t50",
    ]
    assert _convert(tmp_path, lines) == EXPECTED


def test_empty_input_pads_every_chromosome(tmp_path):
    assert _convert(tmp_path, []) == [
        "chr1\t0\t1\t0", "chr10\t0\t1\t0", "chr2\t0\t1\t0", "chr3\t0\t1\t0",
    ]
