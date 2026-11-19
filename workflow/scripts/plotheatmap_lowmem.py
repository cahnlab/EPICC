#!/usr/bin/env python3
"""deeptools plotHeatmap, with its PDF memory kept linear in the column count.

Usage: python3 plotheatmap_lowmem.py <plotHeatmap arguments>

plotHeatmap makes the figure about 4 cm wide per heatmap column and
rasterizes each heatmap (200 dpi by default). In a vector output (PDF, SVG),
matplotlib draws a rasterized artist on an Agg canvas the size of the whole
figure (backend_mixed.MixedModeRenderer.start_rasterizing) and hands the
backend a cropped *view* of it. The PDF backend keeps every image until
savefig writes the file, so each heatmap keeps its full-figure canvas alive:
N canvases, each proportional to N, and peak memory grows with the square of
the number of columns -- about 6 MB x columns^2 at the defaults, past 32 GB
at ~70 columns whatever the region count.

Handing the backend a copy of the crop lets each canvas go as soon as the
heatmap is drawn. The PDF is unchanged: the embedded image streams are byte
for byte the same. With 60 columns of 2,000 regions, peak RSS drops from
22.6 GB to 0.8 GB, and with 120 columns from 88.8 GB to 1.4 GB. What remains
grows with regions x columns: deeptools' own copies of the matrix.
"""

import sys

import numpy as np
from matplotlib.backends import backend_pdf

_draw_image = backend_pdf.RendererPdf.draw_image


def draw_image(self, gc, x, y, im, *args, **kwargs):
    # A copy owns its memory, so the canvas it was cropped from can be freed.
    return _draw_image(self, gc, x, y, np.array(im, copy=True), *args, **kwargs)


backend_pdf.RendererPdf.draw_image = draw_image

if __name__ == "__main__":
    from deeptools.plotHeatmap import main
    sys.exit(main())
