#!/usr/bin/env python3
"""Before/after comparison for plots moved from deeptools (or an older R
script) to the R plotting scripts.

For a plot family, this renders the old and the new version from the same
fixture input and writes:

  <out>/<family>/old.pdf, new.pdf           the two renders
  <out>/<family>/compare.png                old | new side by side, labelled
  <out>/<family>/report.txt                 sizes and whether the renders are
                                            pixel-identical

"old" is the plot as drawn at a git ref (--old-ref: the branch the
conversion PR targets), taken from `git archive` of that ref's workflow/scripts so its
sibling files are the ones it shipped with. "new" is the working tree.
Additional variants (e.g. another theme) are rendered next to "new".

Usage:
  python dev/plot_compare.py --out cmp/ --old-ref REF [--family mc_profile]

Needs Rscript (ggplot2, jsonlite, data.table), pdftoppm (poppler) and Pillow;
the `scales` and `profiles` families also need deeptools on PATH (the old
render). `profiles` runs plotting_profile_on_targetfile's own shell block at
the old ref and now, one compare.png per variant under <out>/profiles/. Each new
conversion phase adds a family to FAMILIES; dev/docs/r-plotting-plan.md lists
which family covers which plot.
"""

import argparse
import hashlib
import json
import math
import subprocess
import sys
import tarfile
import tempfile
from io import BytesIO
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "workflow" / "scripts"


def run(cmd, **kw):
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if proc.returncode != 0:
        sys.exit(f"command failed: {' '.join(map(str, cmd))}\n{proc.stderr}")
    return proc


def scripts_at(ref, dest):
    """Extract workflow/scripts (and the rule files) as of `ref` into dest;
    return the scripts dir."""
    data = subprocess.run(["git", "-C", str(REPO), "archive", ref, "workflow/scripts",
                           "workflow/rules"], capture_output=True, check=True).stdout
    with tarfile.open(fileobj=BytesIO(data)) as tar:
        tar.extractall(dest, filter="data")
    return Path(dest) / "workflow" / "scripts"


def rasterise(pdf, dpi=100):
    prefix = pdf.with_suffix("")
    run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile", str(pdf), str(prefix)])
    return prefix.with_suffix(".png")


def side_by_side(images, out):
    """images: list of (title, png path). Writes one labelled row."""
    from PIL import Image, ImageDraw, ImageFont
    tiles = [(t, Image.open(p).convert("RGB")) for t, p in images]
    pad, head = 20, 40
    width = sum(im.width for _, im in tiles) + pad * (len(tiles) + 1)
    height = max(im.height for _, im in tiles) + head + pad
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 22)
    except OSError:
        font = ImageFont.load_default()
    x = pad
    for title, im in tiles:
        draw.text((x, 8), title, fill="black", font=font)
        canvas.paste(im, (x, head))
        draw.rectangle([x - 1, head - 1, x + im.width, head + im.height], outline="grey")
        x += im.width + pad
    canvas.save(out)


def pixels_md5(png):
    from PIL import Image
    return hashlib.md5(Image.open(png).convert("RGB").tobytes()).hexdigest()


# --- fixtures ------------------------------------------------------------------

def scales_fixture(path, nregions=300, nbins=30):
    """A `most`-style matrix: two marks two orders of magnitude apart, H3K9me
    next to H3K9me2, and a stranded RNA track whose plus/minus columns
    interleave by sample."""
    import gzip
    import random
    rnd = random.Random(7)
    keys = ["WT__H3K9me2", "ddm1__H3K9me2", "WT__H3K9me", "WT__H3K4me3", "ddm1__H3K4me3",
            "WT__RNA_plus", "WT__RNA_minus", "ddm1__RNA_plus", "ddm1__RNA_minus"]
    scale = {"H3K9me2": 60, "H3K9me": 6, "H3K4me3": 0.5, "RNA_plus": 30, "RNA_minus": 3}

    def values(key, r):
        mark = key.split("__")[-1]
        hump = [math.exp(-((b - nbins / 2) / (nbins / 5)) ** 2) for b in range(nbins)]
        level = scale[mark] * (0.5 if key.startswith("ddm1") else 1) * (r / nregions)
        return [max(0.0, level * (0.3 + h) + rnd.gauss(0, scale[mark] * 0.05)) for h in hump]

    half = nregions // 2
    header = {
        "upstream": [1000] * len(keys), "downstream": [1000] * len(keys),
        "body": [1000] * len(keys), "bin size": [100] * len(keys),
        "ref point": [None] * len(keys), "verbose": False, "bin avg type": "mean",
        "missing data as zero": True, "min threshold": None, "max threshold": None,
        "scale": 1, "skip zeros": False, "nan after end": False, "proc number": 1,
        "sort regions": "keep", "sort using": "mean",
        "unscaled 5 prime": [0] * len(keys), "unscaled 3 prime": [0] * len(keys),
        "group_labels": ["genes", "genes_nostrand"],
        "group_boundaries": [0, half, nregions], "sample_labels": keys,
        "sample_boundaries": [i * nbins for i in range(len(keys) + 1)],
    }
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header, separators=(",", ":")) + "\n")
        for r in range(nregions):
            vals = [v for k in keys for v in values(k, r)]
            fh.write(f"chr1\t{r * 5000}\t{r * 5000 + 1000}\tg{r}\t0\t{'+-'[r % 2] if r < half else '.'}\t"
                     + "\t".join(f"{v:.4f}" for v in vals) + "\n")


def family_scales(work, old_scripts, old_smk):
    """heatmaps_scales: type -- computing_matrix_scales as it was (bash +
    deeptools dataRange/plotProfile, frozen in dev/reference/) vs the rule as
    it now stands (R_matrix_scales.R), each drawn by plotHeatmap on the same
    matrix. flags.txt lists both and the largest difference per flag. Needs
    deeptools on PATH."""
    sys.path.insert(0, str(REPO / "dev"))
    from scales_rule import parse_flags, run_scales_rule
    matrix = work / "most.gz"
    scales_fixture(matrix)
    heat, prof = run_scales_rule(matrix, work / "old_rule", "type", "mean")
    smk = (REPO / "workflow" / "rules" / "combined_analysis.smk").read_text()
    new_heat, new_prof = run_scales_rule(matrix, work / "new_rule", "type", "mean", smk_text=smk)
    diffs = []
    for o, n in ((heat, new_heat), (prof, new_prof)):
        o, n = parse_flags(o), parse_flags(n)
        for flag in o:
            d = max(abs(a - b) / max(abs(a), 1e-300) for a, b in zip(o[flag], n[flag]))
            diffs.append(f"  {flag}: largest relative difference {d:.2e}")
    (work / "flags.txt").write_text(
        f"old (bash + dataRange/plotProfile):\n  {heat}\n  {prof}\n"
        f"new (R_matrix_scales.R):\n  {new_heat}\n  {new_prof}\n" + "\n".join(diffs) + "\n")
    renders = []
    for name, h, p in (("old", heat, prof), ("new", new_heat, new_prof)):
        out = work / f"{name}.pdf"
        run(["plotHeatmap", "-m", str(matrix), "-o", str(out), "--colorMap", "seismic",
             "--startLabel", "start", "--endLabel", "end", *h.split(), *p.split()])
        renders.append((f"{name}: " + ("bash + deeptools scales" if name == "old" else "R scales"), out))
    return renders


# --- profiles (phase 2) --------------------------------------------------------

PROFILE_COLOURS = {"WT": "#1f77b4", "ddm1": "#ff7f0e", "WT_EMseq": "#1f77b4"}


def mc_matrix_fixture(path, nregions=240, nbins=30):
    """An mC matrix: sample x context columns with NaNs, a stranded group and
    a nostrand group, two samples sharing a colour (WT, WT_EMseq)."""
    import gzip
    import random
    rnd = random.Random(11)
    series = ["WT", "ddm1", "WT_EMseq"]
    level = {"mCG": 45.0, "mCHG": 20.0, "mCHH": 4.0}
    keys = [f"{s}__{c}" for c in level for s in series]
    half = nregions * 3 // 4
    header = _matrix_header(keys, nbins, ["all_TEs", "all_TEs_nostrand"], [0, half, nregions])
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header, separators=(",", ":")) + "\n")
        for r in range(nregions):
            vals = []
            for key in keys:
                s, c = key.split("__")
                base = level[c] * (0.4 if s == "ddm1" else 1.0) * (0.8 if r >= half else 1.0)
                for b in range(nbins):
                    v = base * (1 + 0.5 * math.exp(-((b - nbins / 2) / (nbins / 5)) ** 2)) \
                        + rnd.gauss(0, base * 0.2) + (0.3 if s == "WT_EMseq" else 0)
                    vals.append("nan" if rnd.random() < 0.05 else f"{max(v, 0):.4f}")
            fh.write(f"chr1\t{r * 5000}\t{r * 5000 + 1000}\tte{r}\t0\t"
                     f"{'+-'[r % 2] if r < half else '.'}\t" + "\t".join(vals) + "\n")
    return keys


def chip_tss_fixture(path, nregions=200, nbins=40):
    """A reference-point (tss) ChIP matrix with one region group."""
    import gzip
    import random
    rnd = random.Random(5)
    keys = ["WT__H3K4me3", "ddm1__H3K4me3", "WT__H3K27me3"]
    header = _matrix_header(keys, nbins, ["all_genes"], [0, nregions], body=0,
                            ref_point="TSS")
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header, separators=(",", ":")) + "\n")
        for r in range(nregions):
            vals = []
            for key in keys:
                lv = {"H3K4me3": 2.0, "H3K27me3": 0.4}[key.split("__")[1]]
                lv *= 0.6 if key.startswith("ddm1") else 1
                for b in range(nbins):
                    vals.append(f"{max(0, lv * math.exp(-((b - 22) / 6) ** 2) + rnd.gauss(0, lv * 0.3)):.4f}")
            fh.write(f"chr1\t{r * 5000}\t{r * 5000 + 1}\tg{r}\t0\t{'+-'[r % 2]}\t"
                     + "\t".join(vals) + "\n")
    return keys


def _matrix_header(keys, nbins, groups, bounds, body=1000, ref_point=None):
    n = len(keys)
    return {
        "upstream": [1000] * n, "downstream": [1000] * n, "body": [body] * n,
        "bin size": [100 if body else 50] * n, "ref point": [ref_point] * n,
        "verbose": False, "bin avg type": "mean", "missing data as zero": False,
        "min threshold": None, "max threshold": None, "scale": 1, "skip zeros": False,
        "nan after end": False, "proc number": 1, "sort regions": "keep", "sort using": "mean",
        "unscaled 5 prime": [0] * n, "unscaled 3 prime": [0] * n,
        "group_labels": groups, "group_boundaries": bounds, "sample_labels": keys,
        "sample_boundaries": [i * nbins for i in range(n + 1)],
    }


def _display(keys):
    rows = []
    for k in keys:
        sample, track = k.split("__")
        mc = track.startswith("mC")
        label = f"{sample}_{track}"
        rows.append(f"{k}\t{label}\t{PROFILE_COLOURS[sample]}\t{track}\t{sample if mc else label}")
    return "\n".join(rows)


def _profile_rule(smk, scripts, work, name, matrix, kind, env, mode, average, display, theme):
    """Run plotting_profile_on_targetfile's shell block from `smk` (the old ref's
    or the working tree's), feeding it the scales files of the same ref's
    computing_matrix_scales, as Snakemake would. Returns the two PDFs."""
    import os
    sys.path.insert(0, str(REPO / "dev"))
    from scales_rule import rule_shell, run_rule_shell, run_scales_rule
    d = work / name
    d.mkdir(parents=True, exist_ok=True)
    run_scales_rule(matrix, d / "scales", mode, average, smk_text=smk, scripts_dir=scripts)
    out = {"plot1": str(d / "Profile.pdf"), "plot2": str(d / "Profile_pergroup.pdf")}
    params = dict(
        analysis_name="cmp", ref_genome="G", target_name="T", matrix=kind, env=env,
        plot_params="--plotType lines", average=average, scales=mode,
        cg_scale=100, chg_scale=50, chh_scale=20, display=display,
        sample_labeller=str(scripts / "matrix_sample_labels.py"),
        mc_script=str(scripts / "R_mc_profile_by_context.R"),
        profile_script=str(scripts / "R_profile.R"), theme=theme)
    inputs = dict(matrix=str(matrix),
                  params_regions=str(d / "scales" / "params_regions.txt"),
                  params_profile=str(d / "scales" / "params_profile.txt"))
    tmp = d / "tmp"
    tmp.mkdir(exist_ok=True)
    run_rule_shell(rule_shell(smk, "plotting_profile_on_targetfile"), params, inputs, out,
                   d / "rule.log", env={**os.environ, "TMPDIR": str(tmp)})
    return Path(out["plot1"]), Path(out["plot2"])


PROFILE_VARIANTS = [
    # (variant, fixture, kind, env, heatmaps_scales, profiles_scale)
    ("most_type_mean", "most", "regions", "most", "type", "mean"),
    ("most_sample_median", "most", "regions", "most", "sample", "median"),
    ("most_default_mean", "most", "regions", "most", "default", "mean"),
    ("mC_type_mean", "mc", "regions", "mC", "type", "mean"),
    ("mC_default_median", "mc", "regions", "mC", "default", "median"),
    ("ChIP_tss_type_mean", "chip", "tss", "ChIP", "type", "mean"),
]


def family_profiles(work, old_scripts, old_smk):
    """Profile__ and Profile_pergroup__: plotting_profile_on_targetfile at the
    old ref (plotProfile; the mC per-group plot from the R dump script) vs
    now (R_profile.R), on the same matrices, for each variant. The old rule
    also gets its scales from the old ref's computing_matrix_scales. Returns
    {variant: renders}; most_type_mean adds the new rule under profile_type
    se and fill and the presentation theme, mC_type_mean under fill and std
    (bands and fill are clipped to each panel's limits)."""
    theme = json.dumps({"overrides": {}, "theme": "default"}, sort_keys=True)
    fixtures = {"most": work / "most.gz", "mc": work / "mc.gz", "chip": work / "chip.gz"}
    scales_fixture(fixtures["most"])
    keys = {"most": ["WT__H3K9me2", "ddm1__H3K9me2", "WT__H3K9me", "WT__H3K4me3",
                     "ddm1__H3K4me3", "WT__RNA_plus", "WT__RNA_minus", "ddm1__RNA_plus",
                     "ddm1__RNA_minus"],
            "mc": mc_matrix_fixture(fixtures["mc"]),
            "chip": chip_tss_fixture(fixtures["chip"])}
    new_smk = (REPO / "workflow" / "rules" / "combined_analysis.smk").read_text()
    out = {}
    for variant, fx, kind, env, mode, average in PROFILE_VARIANTS:
        display = _display(keys[fx])
        old = _profile_rule(old_smk, old_scripts, work, f"{variant}/old", fixtures[fx], kind, env,
                            mode, average, display, theme)
        new = _profile_rule(new_smk, SCRIPTS, work, f"{variant}/new", fixtures[fx], kind, env,
                            mode, average, display, theme)
        extra = []
        styles = {"most_type_mean": (("profile_type se", {"profile_type": "se"}),
                                     ("profile_type fill", {"profile_type": "fill"}),
                                     ("presentation", None)),
                  "mC_type_mean": (("profile_type fill", {"profile_type": "fill"}),
                                   ("profile_type std", {"profile_type": "std"}))}
        if variant in styles:
            for label, spec in styles[variant]:
                t = json.dumps({"overrides": spec or {}, "theme": "presentation" if spec is None
                                else "default"}, sort_keys=True)
                tag = label.replace(" ", "_")
                extra.append((label, _profile_rule(new_smk, SCRIPTS, work, f"{variant}/{tag}",
                                                   fixtures[fx], kind, env, mode, average,
                                                   display, t)))
        out[f"{variant}/Profile"] = [("old (plotProfile)", old[0]), ("new (R)", new[0])] + \
            [(f"new ({lab})", p[0]) for lab, p in extra]
        old_pg = "old (R dump script)" if env == "mC" else "old (plotProfile --perGroup)"
        out[f"{variant}/Profile_pergroup"] = [(old_pg, old[1]), ("new (R)", new[1])] + \
            [(f"new ({lab})", p[1]) for lab, p in extra]
    return out


FAMILIES = {
    "scales": family_scales,
    "profiles": family_profiles,
}


def main():
    ap = argparse.ArgumentParser(description="Before/after comparison of converted plots")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--family", choices=sorted(FAMILIES), action="append")
    ap.add_argument("--old-ref", required=True,
                    help="git ref for the old render, normally the PR's base branch")
    a = ap.parse_args()
    ref = a.old_ref
    with tempfile.TemporaryDirectory() as tmp:
        old_scripts = scripts_at(ref, tmp)
        for fam in a.family or sorted(FAMILIES):
            work = a.out / fam
            work.mkdir(parents=True, exist_ok=True)
            old_smk = (old_scripts.parent / "rules" / "combined_analysis.smk").read_text()
            result = FAMILIES[fam](work, old_scripts, old_smk)
            # A family returns one list of renders, or {variant: renders}.
            variants = result if isinstance(result, dict) else {"": result}
            report = [f"family: {fam}", f"old ref: {ref}"]
            for variant, renders in variants.items():
                vdir = work / variant if variant else work
                vdir.mkdir(parents=True, exist_ok=True)
                pngs = [(t, rasterise(p)) for t, p in renders]
                side_by_side(pngs, vdir / "compare.png")
                old_md5 = pixels_md5(pngs[0][1])
                if variant:
                    report.append(f"[{variant}] -> {vdir / 'compare.png'}")
                for (t, p), (_, png) in zip(renders, pngs):
                    same = "identical to old" if pixels_md5(png) == old_md5 else "differs from old"
                    report.append(f"  {t}: {p.name}, {p.stat().st_size} bytes, pixels {same}")
            (work / "report.txt").write_text("\n".join(report) + "\n")
            print("\n".join(report), "\n")


if __name__ == "__main__":
    main()
