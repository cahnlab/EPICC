"""Run computing_matrix_scales' shell block on a matrix, old or new.

  run_scales_rule(...)            the old bash + deeptools version, frozen in
                                  dev/reference/computing_matrix_scales_deeptools.sh
                                  (needs computeMatrixOperations and plotProfile)
  run_scales_rule(..., smk_text=) the block in a .smk text (the working tree's:
                                  Rscript R_matrix_scales.R)

The block's Snakemake placeholders are filled in the way Snakemake fills them
(str.format with params/input/output/log namespaces) and it is run with bash.
"""

import re
import shlex
import string
import subprocess
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]


REFERENCE = REPO / "dev" / "reference" / "computing_matrix_scales_deeptools.sh"


def reference_shell():
    return REFERENCE.read_text()


def scales_shell(smk_text):
    start = smk_text.index("rule computing_matrix_scales:")
    m = re.compile(r'shell:\s*"""(.*?)"""', re.S).search(smk_text, start)
    return m.group(1)


class _SnakemakeFormatter(string.Formatter):
    """str.format plus Snakemake's ':q' (shell-quote) format spec."""

    def format_field(self, value, spec):
        if spec == "q":
            return shlex.quote(str(value))
        return super().format_field(value, spec)


def rule_shell(smk_text, rule):
    """The shell block of `rule` in a .smk text."""
    start = smk_text.index(f"rule {rule}:")
    m = re.compile(r'shell:\s*"""(.*?)"""', re.S).search(smk_text, start)
    return m.group(1)


def run_rule_shell(shell, params, inputs, outputs, log, env=None):
    """Fill a rule's shell block the way Snakemake does and run it with bash
    (TMPDIR as given in env, or the caller's)."""
    script = _SnakemakeFormatter().format(
        shell, params=SimpleNamespace(**params), input=SimpleNamespace(**inputs),
        output=SimpleNamespace(**outputs), log=str(log), wildcards=SimpleNamespace())
    proc = subprocess.run(["bash", "-c", "set -euo pipefail\n" + script],
                          capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"rule failed:\n{proc.stdout}\n{proc.stderr}")
    return proc


def run_scales_rule(matrix, workdir, scales="type", average="mean",
                    mc=(100, 50, 20), smk_text=None, scripts_dir=None):
    """Returns (params_heatmap text, params_profile text). smk_text=None runs
    the frozen deeptools reference; otherwise the block in that .smk text."""
    shell = reference_shell() if smk_text is None else scales_shell(smk_text)
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
        scales_script=str(scripts_dir / "R_matrix_scales.R"),
        cg_scale=mc[0], chg_scale=mc[1], chh_scale=mc[2])
    script = shell.format(
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
