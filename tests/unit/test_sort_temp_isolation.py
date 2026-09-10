"""samtools sort temp files must live in the per-job $TMPDIR.

samtools sort refuses to overwrite an existing temp chunk, drops the reads it
could not write, and still exits 0. A prefix inside the results tree therefore
turns any killed attempt into silent read loss on the retry, with no failing
exit status anywhere to catch it.
"""

import glob
import os
import re
import shutil
import subprocess

import pytest

_RULES = sorted(glob.glob(os.path.join(
    os.path.dirname(__file__), "..", "..", "workflow", "rules", "*.smk")))

# A -T prefix is safe if it sits under $TMPDIR or a mktemp -d dir.
_SAFE = re.compile(r'-T\s+"?\$(TMPDIR|tmpd)')


def _sort_calls(text):
    for line in text.splitlines():
        if "samtools sort" in line and "-T" in line:
            yield line.strip()


@pytest.mark.parametrize("path", _RULES, ids=lambda p: os.path.basename(p))
class TestSortPrefixes:
    def test_every_prefix_is_job_local(self, path):
        with open(path) as fh:
            offenders = [c for c in _sort_calls(fh.read())
                         if not _SAFE.search(c)]
        assert not offenders, (
            f"{os.path.basename(path)}: sort temp prefix outside $TMPDIR:\n  "
            + "\n  ".join(offenders))

    def test_no_prefix_derives_from_an_output_path(self, path):
        # -T {output.x}.sort puts the chunks next to the output, where a
        # killed attempt leaves them for the retry to trip over.
        with open(path) as fh:
            offenders = [c for c in _sort_calls(fh.read())
                         if re.search(r"-T\s+\"?\{output", c)]
        assert not offenders, "\n  ".join(offenders)


class TestChipFilterIntermediate:
    """The pre-markdup BAM is untracked, so snakemake cannot clean it up."""

    def test_not_written_into_the_results_tree(self):
        path = os.path.join(os.path.dirname(__file__), "..", "..",
                            "workflow", "rules", "ChIPseq.smk")
        with open(path) as fh:
            text = fh.read()
        assert "sorted_{params.sample_name}" not in text
        assert text.count('sorted_bam="$TMPDIR/sorted.bam"') == 2


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools missing")
class TestSamtoolsSortBehaviour:
    """Why the prefix has to be clean: sort lies about failing."""

    @staticmethod
    def _make_bam(tmp_path, n=200000):
        sam = tmp_path / "in.sam"
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:20000000\n")
            for i in range(1, n + 1):
                fh.write(f"r{i}\t0\tchr1\t{(n - i) * 20 + 1}\t60\t10M\t*\t0\t0"
                         "\tACGTACGTAC\tIIIIIIIIII\n")
        bam = tmp_path / "in.bam"
        subprocess.run(["samtools", "view", "-bo", str(bam), str(sam)],
                       check=True, capture_output=True)
        return bam, n

    @staticmethod
    def _count(bam):
        out = subprocess.run(["samtools", "view", "-c", str(bam)],
                             check=True, capture_output=True, text=True)
        return int(out.stdout.strip())

    def test_a_stale_chunk_costs_reads_without_failing(self, tmp_path):
        bam, n = self._make_bam(tmp_path)
        stale = tmp_path / "p.sort.0000.bam"
        subprocess.run(["samtools", "view", "-bo", str(stale), str(bam)],
                       check=True, capture_output=True)
        out = tmp_path / "out.bam"
        # -m 1M is the documented floor and forces chunks to spill to disk.
        result = subprocess.run(
            ["samtools", "sort", "-m", "1M", "-T", str(tmp_path / "p.sort"),
             "-o", str(out), str(bam)], capture_output=True, text=True)
        assert result.returncode == 0, "guard is moot if sort ever fails loudly"
        assert "File exists" in result.stderr
        assert self._count(out) < n, "expected silent read loss"

    def test_a_clean_prefix_keeps_every_read(self, tmp_path):
        bam, n = self._make_bam(tmp_path)
        out = tmp_path / "out.bam"
        subprocess.run(
            ["samtools", "sort", "-m", "1M", "-T", str(tmp_path / "clean.sort"),
             "-o", str(out), str(bam)], check=True, capture_output=True)
        assert self._count(out) == n
