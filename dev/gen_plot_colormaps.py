#!/usr/bin/env python3
"""Write workflow/scripts/plot_colormaps.json: matplotlib's colour maps as hex
stops, so a theme (or a translated `--colorMap NAME` from the deprecated
heatmaps_plot_params) can name a colour map the way deeptools took one.

plot_theme.R resolves a one-element `heatmap_colors` entry through this table
(a `_r` suffix reverses the map); config_validation.py checks names against it.
17 evenly spaced stops keep the quarter points of the diverging maps (seismic,
bwr, ...) exact and are plenty for a continuous ramp.

Run with any Python that has matplotlib (the epibutton env does):
  python dev/gen_plot_colormaps.py
"""

import json
from pathlib import Path

import matplotlib
from matplotlib.colors import to_hex

STOPS = 17
OUT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "plot_colormaps.json"


def main():
    maps = {}
    for name in sorted(n for n in matplotlib.colormaps if not n.endswith("_r")):
        cmap = matplotlib.colormaps[name]
        maps[name] = [to_hex(cmap(i / (STOPS - 1))) for i in range(STOPS)]
    source = (f"matplotlib {matplotlib.__version__}, {STOPS} evenly spaced stops; "
              "regenerate with dev/gen_plot_colormaps.py")
    # One map per line, so a regenerated table diffs by map.
    lines = [f"  {json.dumps(n)}: {json.dumps(maps[n])}" for n in maps]
    OUT.write_text("{\n" + f'"_source": {json.dumps(source)},\n"maps": {{\n'
                   + ",\n".join(lines) + "\n}\n}\n")
    print(f"{len(maps)} colour maps -> {OUT}")


if __name__ == "__main__":
    main()
