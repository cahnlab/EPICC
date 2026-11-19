"""Memory request for the plotHeatmap rules, from the matrix they plot.

plotHeatmap runs through plotheatmap_lowmem.py, which keeps the PDF's
rasterized heatmaps from growing with the square of the column count. What is
left grows with the matrix: deeptools reads it as masked float64 rows,
stacks them and copies them again to sort and clip, about 31 bytes per value
in total; each heatmap column also keeps its rasterized panel (~3 MB at the
default size and 200 dpi) until the PDF is written. Measured peak RSS
(deeptools 4.0.0, 180 bins per column):

    regions x columns   values   peak
    2,000 x 60          22 M     0.8 GB
    2,000 x 120         43 M     1.4 GB
    20,000 x 30         108 M    3.5 GB
    50,000 x 60         540 M    16.8 GB

The request is that model with a 1.5x margin, times the attempt number so a
retry asks for more. Snakemake evaluates it again just before it schedules
the job, when the matrix exists; until then (DAG building, dry-runs) the
matrix is missing and the old fixed request is used.
"""

import gzip
import json
import math

BYTES_PER_VALUE = 32
BASE_MB = 1000
PANEL_MB = 3
MARGIN = 1.5
MIN_MB = 4000
FALLBACK_MB = 32000


def matrix_dimensions(matrix_path):
    """(regions, total bins over all columns, columns) from the matrix header."""
    with gzip.open(matrix_path, "rt") as fh:
        header = json.loads(fh.readline().lstrip("@"))
    return (header["group_boundaries"][-1], header["sample_boundaries"][-1],
            len(header["sample_labels"]))


def heatmap_mem_mb(matrix_path, attempt=1):
    """mem_mb for plotting the heatmap of `matrix_path`, rounded up to 1000."""
    try:
        regions, bins, columns = matrix_dimensions(matrix_path)
    except (OSError, ValueError, KeyError, IndexError, EOFError):
        return FALLBACK_MB * attempt
    mb = BASE_MB + regions * bins * BYTES_PER_VALUE / 1e6 + columns * PANEL_MB
    mb = max(MIN_MB, mb * MARGIN) * attempt
    return int(math.ceil(mb / 1000) * 1000)
