"""Run the bash of combined_analysis.smk's computing_matrix_scales on a matrix.

Used to compare plot_scales.R with the rule it replaces: the shell block is
taken from the .smk text (the working tree's or a git ref's), its Snakemake
placeholders are filled in the way Snakemake fills them (str.format with
params/input/output/log namespaces), and it is run with bash. Needs deeptools
(computeMatrixOperations, plotProfile) on PATH.
"""

import re
import subprocess
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]


def scales_shell(smk_text):
    start = smk_text.index("rule computing_matrix_scales:")
    m = re.compile(r'shell:\s*"""(.*?)"""', re.S).search(smk_text, start)
    return m.group(1)


def run_scales_rule(matrix, workdir, scales="type", average="mean",
                    mc=(100, 50, 20), smk_text=None, scripts_dir=None):
    """Returns (params_heatmap text, params_profile text)."""
    if smk_text is None:
        smk_text = (REPO / "workflow" / "rules" / "combined_analysis.smk").read_text()
    scripts_dir = Path(scripts_dir or REPO / "workflow" / "scripts")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out = {k: str(workdir / f"{k}.txt") for k in
           ("params_heatmap", "params_profile", "params_regions", "temp_values",
            "temp_profile_values")}
    out["temp_profile"] = str(workdir / "temp_profile.pdf")
    params = SimpleNamespace(
        analysis_name="cmp", ref_genome="G", env="most", target_name="T", matrix="regions",
        scales=scales, profile=average,
        labeller=str(scripts_dir / "matrix_region_labels.py"),
        cg_scale=mc[0], chg_scale=mc[1], chh_scale=mc[2])
    script = scales_shell(smk_text).format(
        params=params, input=SimpleNamespace(matrix=str(matrix)),
        output=SimpleNamespace(**out), log=str(workdir / "rule.log"),
        wildcards=SimpleNamespace())
    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"scales rule failed:\n{proc.stdout}\n{proc.stderr}")
    return (Path(out["params_heatmap"]).read_text(), Path(out["params_profile"]).read_text())


def parse_flags(text):
    """'--zMin 1 2 --zMax 3 4' -> {'--zMin': [1.0, 2.0], '--zMax': [3.0, 4.0]}."""
    flags, cur = {}, None
    for w in text.split():
        if w.startswith("--"):
            cur = flags.setdefault(w, [])
        else:
            cur.append(float(w))
    return flags
