"""The stranded halves of a deeptools matrix must merge into ONE region group.

`making_stranded_matrix_on_targetfile` splits a stranded target file into a plus
and a minus BED and runs computeMatrix on each; `merging_matrix` rbinds the two
back together. Every rule downstream of that -- plotHeatmap, plotProfile, and
sort_heatmap's `relabel --groupLabels` -- passes exactly one `--regionsLabel`,
so the merged matrix has to carry exactly one group.

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
from pathlib import Path

import pytest

pytestmark = pytest.mark.requires_deeptools

if shutil.which("computeMatrixOperations") is None:
    pytest.skip("deeptools not found", allow_module_level=True)


REPO = Path(__file__).resolve().parents[2]
SMK = REPO / "workflow" / "rules" / "combined_analysis.smk"

N_PLUS = 6
N_MINUS = 4
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

    def test_unstranded_input_passes_through(self, tmp_path):
        """One input means an unstranded target file: copied, still one group."""
        only = tmp_path / "matrix_tss__unstranded.gz"
        _write_matrix(only, "temp_file_tss__mC__all_genes_unstranded", N_PLUS)
        merged = _run_merge([only], tmp_path)
        header = _read_header(merged)
        assert len(header["group_labels"]) == 1
        assert header["group_boundaries"] == [0, N_PLUS]
