"""Read the JSON header of a computeMatrix matrix file.

A deeptools matrix is a gzipped file whose first line is '@' followed by a JSON
header (group and sample labels and boundaries, bin layout); the region rows
follow. Shared by matrix_region_labels.py and matrix_sample_labels.py, which
run as standalone scripts, so this sits beside them and is imported by name.
"""

import gzip
import json


def read_matrix_header(matrix_path):
    with gzip.open(matrix_path, "rt") as fh:
        return json.loads(fh.readline().lstrip("@"))
