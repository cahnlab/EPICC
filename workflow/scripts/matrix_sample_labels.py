#!/usr/bin/env python3
"""Print the display row of every sample column in a computeMatrix matrix.

Matrices name their sample columns by stable keys (define_key_for_plots), not
by plot labels, so that a label change reruns only the plotting jobs and never
rebuilds a matrix. The plotting rules pass their labels and colours as a
display table keyed the same way; this looks each column of the matrix up in
it, in column order, ready for plotHeatmap/plotProfile --samplesLabel and
--colors.

Reading the order from the matrix rather than assuming it means the stranded
(plus/minus/nostrand) and unstranded matrices, which carry different columns,
both come out right from one display table.

Usage: matrix_sample_labels.py <matrix.gz> <display.tsv>
  display.tsv: key<TAB>label<TAB>... per line (extra columns passed through)
  -> the display line of each matrix column, in column order
"""

import gzip
import json
import sys


def sample_keys(matrix_path):
    with gzip.open(matrix_path, "rt") as fh:
        header = json.loads(fh.readline().lstrip("@"))
    return header["sample_labels"]


def read_display(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                rows[line.split("\t", 1)[0]] = line
    return rows


def display_rows(matrix_path, display_path):
    rows = read_display(display_path)
    keys = sample_keys(matrix_path)
    missing = [k for k in keys if k not in rows]
    if missing:
        raise SystemExit(
            f"{matrix_path}: no label for sample column(s) {', '.join(missing)}. "
            "The matrix predates the current sample sheet; delete it to rebuild."
        )
    return [rows[k] for k in keys]


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    for row in display_rows(sys.argv[1], sys.argv[2]):
        print(row)
