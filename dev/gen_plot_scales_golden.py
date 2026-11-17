#!/usr/bin/env python3
"""Write tests/unit/data/plot_scales_golden.json: what the old
computing_matrix_scales (bash + deeptools dataRange/plotProfile, frozen in
dev/reference/computing_matrix_scales_deeptools.sh) wrote for the fixture
matrix of tests/unit/test_plot_scales.py, for every heatmaps_scales mode and
profiles_scale average. test_plot_scales.py checks the rule as it now stands
against it without needing deeptools.

Needs deeptools on PATH. Run from the repository root:
  python dev/gen_plot_scales_golden.py
"""

import gzip
import hashlib
import json
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "dev"))
from scales_rule import run_scales_rule  # noqa: E402
from tests.unit.test_plot_scales import MC, GOLDEN, build  # noqa: E402


def main():
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        matrix = Path(tmp) / "m.gz"
        build(matrix)
        with gzip.open(matrix, "rb") as fh:
            out["fixture_sha256"] = hashlib.sha256(fh.read()).hexdigest()
        for mode in ("type", "sample", "default"):
            for average in ("mean", "median"):
                out[f"{mode}_{average}"] = list(run_scales_rule(
                    matrix, Path(tmp) / f"{mode}_{average}", mode, average, MC))
    GOLDEN.parent.mkdir(parents=True, exist_ok=True)
    GOLDEN.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(f"-> {GOLDEN}")


if __name__ == "__main__":
    main()
