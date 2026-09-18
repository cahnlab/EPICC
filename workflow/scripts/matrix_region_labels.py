#!/usr/bin/env python3
"""Print the --regionsLabel argument for a computeMatrix matrix.

plotHeatmap, plotProfile and computeMatrixOperations relabel all want exactly
one label per region group, and fail outright on a miscount ("length new labels
!= length original labels"). The group structure is recorded in the matrix
itself, so read it from there rather than guessing from the target file: a
stranded target merges its plus and minus halves into one group but keeps a
'nostrand' group of its own, and computeMatrix's --skipZeros drops regions, so
neither the group count nor the region counts can be derived from the BED.

Usage: matrix_region_labels.py <matrix.gz>
  -> --regionsLabel all_TEs(42296) all_TEs_nostrand(631)
"""

import gzip
import json
import sys


def region_labels(matrix_path):
    with gzip.open(matrix_path, "rt") as fh:
        header = json.loads(fh.readline().lstrip("@"))
    bounds = header["group_boundaries"]
    labels = header["group_labels"]
    if len(bounds) != len(labels) + 1:
        raise ValueError(
            f"{matrix_path}: {len(labels)} group labels but {len(bounds)} "
            "boundaries; the matrix header is inconsistent"
        )
    return [f"{lab}({bounds[i + 1] - bounds[i]})" for i, lab in enumerate(labels)]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    # No trailing newline: the caller interpolates this straight into a command.
    sys.stdout.write("--regionsLabel " + " ".join(region_labels(sys.argv[1])))
