"""Memory request for the heatmap rule, from the matrix it plots.

The heatmaps are drawn in R (R_heatmap.R on plot_heatmap.R). The script reads
the matrix once into a data.table (8 bytes per value, plus the decompressed
text fread maps while it parses) and then works one sample column at a time:
each column's rows are averaged down to the pixel rows they are drawn at, so
nothing that scales with the region count is copied a second time, and the
figure itself costs a couple of MB per column. heatmaps_scales: default finds
its percentiles column by column too. Measured peak RSS (R 4.5, ComplexHeatmap
2.x, 180 bins per column, theme defaults; heatmaps_scales type, mean sort):

    regions x columns   values   peak (R)   deeptools 4 + lowmem wrapper
    2,000 x 60          22 M     0.9 GB     0.8 GB
    2,000 x 120         43 M     1.4 GB     1.4 GB
    20,000 x 30         108 M    2.6 GB     3.5 GB
    50,000 x 60         540 M    10.5 GB    16.8 GB

heatmaps_scales: default and heatmaps_sort_options: median measured the same
within 3%. The model is 0.5 GB + 19 bytes per value + 2 MB per column, with a
1.5x margin, at least 4 GB, times the attempt number so a retry asks for more.
A larger heatmap_raster_dpi or heatmap_height in the theme adds to the
per-column part (it grows with the pixels drawn). Snakemake evaluates the
request again just before it schedules the job, when the matrix exists; until
then (DAG building, dry-runs) the matrix is missing and a fixed fallback is
used.
"""

import gzip
import json
import math

BYTES_PER_VALUE = 19
BASE_MB = 500
PANEL_MB = 2
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
