"""merging_matrix must produce exactly the region groups the plots expect.

`making_stranded_matrix_on_targetfile` splits a stranded target file by strand
and runs computeMatrix on each part; `merging_matrix` rbinds them back. The plus
and minus halves are the same regions seen from their own sense strand, so they
have to come back as ONE group. Rows whose strand is '.' or '?' have no sense at
all, so they run as a third pass and stay a group of their own rather than being
folded in or dropped. Everything downstream -- plotHeatmap, plotProfile, and
sort_heatmap's `relabel --groupLabels` -- passes one `--regionsLabel` per group
and fails on a miscount, so the group structure here is load-bearing.

It did, until deeptools 4. deeptools 3 labelled a single-BED region group
'genes' no matter what the file was called, so both halves shared a label and
rbind folded them together for free. deeptools 4 labels the group after the BED
basename, the halves stopped matching, and every heatmap and profile in the
pipeline died on "length new labels != length original labels".

These tests run the real shell body out of combined_analysis.smk rather than a
copy of it, so a future edit that drops the relabel fails here.
"""

import gzip
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.requires_deeptools

if shutil.which("computeMatrixOperations") is None:
    pytest.skip("deeptools not found", allow_module_level=True)


REPO = Path(__file__).resolve().parents[2]
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"

N_PLUS = 6
N_MINUS = 4
N_NOSTRAND = 3
N_BINS = 5


# ---------------------------------------------------------------------------
# Building a minimal deeptools matrix
# ---------------------------------------------------------------------------

def _write_matrix(path, group_label, n_regions, start=0):
    """Write the smallest matrix file deeptools will read back.

    A matrix is a gzipped '@'-prefixed JSON header followed by one BED6+values
    row per region, so it can be built without a bigwig.
    """
    header = {
        "upstream": [200], "downstream": [200], "body": [0], "bin size": [100],
        "ref point": ["TSS"], "verbose": False, "bin avg type": "mean",
        "missing data as zero": False, "min threshold": None,
        "max threshold": None, "scale": 1, "skip zeros": False,
        "nan after end": False, "proc number": 1, "sort regions": "keep",
        "sort using": "mean", "unscaled 5 prime": [0], "unscaled 3 prime": [0],
        "group_labels": [group_label],
        "group_boundaries": [0, n_regions],
        "sample_labels": ["sample1"],
        "sample_boundaries": [0, N_BINS],
    }
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header) + "\n")
        for i in range(n_regions):
            pos = (start + i) * 1000
            row = ["chr1", str(pos), str(pos + 500), f"gene{start + i}", "0", "+"]
            row += [f"{(i + j) / 10:.2f}" for j in range(N_BINS)]
            fh.write("\t".join(row) + "\n")


def _read_header(path):
    with gzip.open(path, "rt") as fh:
        return json.loads(fh.readline().lstrip("@"))


# ---------------------------------------------------------------------------
# Running the rule's own shell body
# ---------------------------------------------------------------------------

def _render_merging_matrix(inputs, output, log, target_name="all_genes"):
    """Extract merging_matrix's shell body and expand it the way Snakemake does.

    Snakemake substitutes the {...} placeholders and then collapses the doubled
    braces that protect literal bash expansions, so do it in that order here.
    """
    text = SMK.read_text()
    start = text.index("rule merging_matrix:")
    block = text[start:text.index("rule computing_matrix_scales:", start)]
    body = re.search(r'shell:\s*"""(.*?)"""', block, re.DOTALL).group(1)

    for placeholder, value in (
        ("{input}", " ".join(str(p) for p in inputs)),
        ("{output}", str(output)),
        ("{log}", str(log)),
        ("{params.matrix}", "tss"),
        ("{params.env}", "mC"),
        ("{params.target_name}", target_name),
        ("{params.ref_genome}", "TestGenome"),
    ):
        body = body.replace(placeholder, value)
    return body.replace("{{", "{").replace("}}", "}")


def _run_merge(inputs, tmp_path, target_name="all_genes"):
    out = tmp_path / "final_matrix.gz"
    log = tmp_path / "merge.log"
    script = _render_merging_matrix(inputs, out, log, target_name)
    workdir = tmp_path / "tmpdir"
    workdir.mkdir()
    proc = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True,
        env={"PATH": __import__("os").environ["PATH"], "TMPDIR": str(workdir)},
    )
    assert proc.returncode == 0, f"merging_matrix failed:\n{proc.stdout}\n{proc.stderr}"
    return out


@pytest.fixture
def stranded_inputs(tmp_path):
    """The two matrices merging_matrix receives for a stranded target file.

    deeptools 4 names each group after its input BED, so the labels differ --
    which is exactly the condition that broke the merge.
    """
    plus = tmp_path / "matrix_tss__plus.gz"
    minus = tmp_path / "matrix_tss__minus.gz"
    _write_matrix(plus, "temp_file_tss__mC__all_genes_plus", N_PLUS, start=0)
    _write_matrix(minus, "temp_file_tss__mC__all_genes_minus", N_MINUS, start=N_PLUS)
    return plus, minus


class TestStrandedMatrixMerge:
    def test_halves_merge_into_one_group(self, stranded_inputs, tmp_path):
        """The regression: two differently-labelled halves, one group out."""
        merged = _run_merge(stranded_inputs, tmp_path)
        header = _read_header(merged)
        assert len(header["group_labels"]) == 1, (
            f"merged matrix has {len(header['group_labels'])} region groups "
            f"({header['group_labels']}); every downstream rule passes a single "
            "--regionsLabel and will fail on the count mismatch"
        )

    def test_group_takes_the_target_name(self, stranded_inputs, tmp_path):
        merged = _run_merge(stranded_inputs, tmp_path, target_name="all_genes")
        assert _read_header(merged)["group_labels"] == ["all_genes"]

    def test_no_regions_are_lost(self, stranded_inputs, tmp_path):
        merged = _run_merge(stranded_inputs, tmp_path)
        assert _read_header(merged)["group_boundaries"] == [0, N_PLUS + N_MINUS]

    def test_single_regions_label_is_accepted(self, stranded_inputs, tmp_path):
        """What plotHeatmap, plotProfile and sort_heatmap all do to the result."""
        merged = _run_merge(stranded_inputs, tmp_path)
        proc = subprocess.run(
            ["computeMatrixOperations", "relabel", "-m", str(merged),
             "--groupLabels", f"all_genes({N_PLUS + N_MINUS})",
             "-o", str(tmp_path / "relabelled.gz")],
            capture_output=True, text=True,
        )
        assert proc.returncode == 0, (
            "a single group label was rejected by the merged matrix:\n"
            f"{proc.stdout}\n{proc.stderr}"
        )

    def test_unstranded_input_is_relabelled_not_copied(self, tmp_path):
        """One input means an unstranded target file: one group, named usefully.

        computeMatrix names the group after its input BED, which on this path is
        an internal temp file -- and that name is what gets printed on the plot.
        """
        only = tmp_path / "matrix_tss__unstranded.gz"
        _write_matrix(only, "temp_file_tss__mC__all_genes_unstranded", N_PLUS)
        merged = _run_merge([only], tmp_path)
        header = _read_header(merged)
        assert header["group_labels"] == ["all_genes"]
        assert header["group_boundaries"] == [0, N_PLUS]


@pytest.fixture
def mixed_strand_inputs(tmp_path):
    """What a target file carrying '.' or '?' rows produces: three passes.

    An annotation with a handful of unstranded rows -- 631 of ColCEN's 42,927
    TEs -- must still be treated as stranded, with those rows held apart rather
    than dropped or folded in, since sense is undefined for them.
    """
    plus = tmp_path / "matrix_tss__plus.gz"
    minus = tmp_path / "matrix_tss__minus.gz"
    nostrand = tmp_path / "matrix_tss__nostrand.gz"
    _write_matrix(plus, "temp_file_tss__mC__all_TEs_plus", N_PLUS, start=0)
    _write_matrix(minus, "temp_file_tss__mC__all_TEs_minus", N_MINUS, start=N_PLUS)
    _write_matrix(nostrand, "temp_file_tss__mC__all_TEs_nostrand", N_NOSTRAND,
                  start=N_PLUS + N_MINUS)
    return plus, minus, nostrand


class TestMixedStrandMerge:
    def test_two_groups_stranded_and_nostrand(self, mixed_strand_inputs, tmp_path):
        merged = _run_merge(mixed_strand_inputs, tmp_path, target_name="all_TEs")
        assert _read_header(merged)["group_labels"] == ["all_TEs", "all_TEs_nostrand"]

    def test_plus_and_minus_still_fold_together(self, mixed_strand_inputs, tmp_path):
        """The nostrand pass must not split plus from minus as a side effect."""
        merged = _run_merge(mixed_strand_inputs, tmp_path, target_name="all_TEs")
        bounds = _read_header(merged)["group_boundaries"]
        assert bounds == [0, N_PLUS + N_MINUS, N_PLUS + N_MINUS + N_NOSTRAND]

    def test_no_regions_are_dropped(self, mixed_strand_inputs, tmp_path):
        merged = _run_merge(mixed_strand_inputs, tmp_path, target_name="all_TEs")
        total = N_PLUS + N_MINUS + N_NOSTRAND
        assert _read_header(merged)["group_boundaries"][-1] == total

    def test_region_labels_match_the_group_count(self, mixed_strand_inputs, tmp_path):
        """What computing_matrix_scales emits must fit what the matrix holds."""
        merged = _run_merge(mixed_strand_inputs, tmp_path, target_name="all_TEs")
        out = subprocess.run(
            [sys.executable, str(REPO / "workflow" / "scripts" / "matrix_region_labels.py"),
             str(merged)], capture_output=True, text=True, check=True).stdout
        assert out == f"--regionsLabel all_TEs({N_PLUS + N_MINUS}) all_TEs_nostrand({N_NOSTRAND})"
        # and the consumers accept it
        labels = out.split()[1:]
        proc = subprocess.run(
            ["computeMatrixOperations", "relabel", "-m", str(merged),
             "--groupLabels", *labels, "-o", str(tmp_path / "relabelled.gz")],
            capture_output=True, text=True,
        )
        assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"

    def test_sort_heatmap_round_trip_keeps_both_groups(self, mixed_strand_inputs, tmp_path):
        """sort_heatmap relabels, then sorts against plotHeatmap's region dump.

        `computeMatrixOperations sort` matches groups by NAME against the
        deepTools_group column of that dump -- which is why the relabel has to
        come first, and why it has to produce every group, not just the first.
        """
        merged = _run_merge(mixed_strand_inputs, tmp_path, target_name="all_TEs")
        labels = subprocess.run(
            [sys.executable, str(REPO / "workflow" / "scripts" / "matrix_region_labels.py"),
             str(merged)], capture_output=True, text=True, check=True).stdout.split()[1:]

        relabelled = tmp_path / "relabelled.gz"
        subprocess.run(
            ["computeMatrixOperations", "relabel", "-m", str(merged),
             "--groupLabels", *labels, "-o", str(relabelled)], check=True,
            capture_output=True,
        )

        # what plotHeatmap --outFileSortedRegions writes
        header = _read_header(relabelled)
        bounds = header["group_boundaries"]
        with gzip.open(relabelled, "rt") as fh:
            rows = [l.rstrip("\n").split("\t")[:6] for l in fh.readlines()[1:]]
        dump = tmp_path / "sorted_regions.bed"
        with open(dump, "w") as fh:
            fh.write("#chrom\tstart\tend\tname\tscore\tstrand\tdeepTools_group\n")
            for i, lab in enumerate(labels):
                for row in rows[bounds[i]:bounds[i + 1]]:
                    fh.write("\t".join(row + [lab]) + "\n")

        out = tmp_path / "sorted.gz"
        proc = subprocess.run(
            ["computeMatrixOperations", "sort", "-m", str(relabelled),
             "-R", str(dump), "-o", str(out)], capture_output=True, text=True,
        )
        assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"
        assert _read_header(out)["group_boundaries"] == [
            0, N_PLUS + N_MINUS, N_PLUS + N_MINUS + N_NOSTRAND
        ]
