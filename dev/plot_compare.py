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
the `scales` family also needs deeptools on PATH. Each new
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

def mc_profile_fixture(d):
    """A plotProfile --outFileNameData dump of an mC matrix: three samples
    (two sharing a colour, so dash patterns show), three contexts, a stranded
    target with a nostrand group."""
    series = ["WT", "ddm1", "WT_EMseq"]
    colours = {"WT": "#1f77b4", "ddm1": "#ff7f0e", "WT_EMseq": "#1f77b4"}
    contexts = ["mCG", "mCHG", "mCHH"]
    level = {"mCG": 45.0, "mCHG": 20.0, "mCHH": 4.0}
    groups = ["all_TEs", "all_TEs_nostrand"]
    nbins = 90
    ticks = [""] * nbins
    ticks[0], ticks[20], ticks[69], ticks[89] = "-2.0Kb", "start", "end", "2.0Kb"
    labels = [f"{s}_{c}" for c in contexts for s in series]
    lines = ["bin labels\t\t" + "\t".join(ticks * len(labels)),
             "bins\t\t" + "\t".join(str(float(i + 1)) for i in range(nbins))]
    for c in contexts:
        for k, s in enumerate(series):
            for gi, g in enumerate(groups):
                base = level[c] * (0.4 if s == "ddm1" else 1.0) * (0.8 if gi else 1.0)
                vals = [base * (1 + 0.5 * math.exp(-((i - 45) / 18) ** 2)) + k * 0.3
                        for i in range(nbins)]
                lines.append(f"{s}_{c}\t{g}\t" + "\t".join(f"{v:.4f}" for v in vals))
    (d / "profile.tab").write_text("\n".join(lines) + "\n")
    (d / "samples.tsv").write_text("".join(
        f"{s}_{c}\t{c}\t{colours[s]}\t{s}\n" for c in contexts for s in series))


def family_mc_profile(work, old_scripts, old_smk):
    """Profile_pergroup__*__mC: R_mc_profile_by_context.R, old ref vs now."""
    mc_profile_fixture(work)
    args = [str(work / "profile.tab"), str(work / "samples.tsv")]
    renders = []
    run(["Rscript", str(old_scripts / "R_mc_profile_by_context.R"), *args,
         str(work / "old.pdf"), "mean methylation"])
    renders.append(("old", work / "old.pdf"))
    run(["Rscript", str(SCRIPTS / "R_mc_profile_by_context.R"), *args,
         str(work / "new.pdf"), "mean methylation"])
    renders.append(("new (default theme)", work / "new.pdf"))
    theme = work / "presentation.json"
    theme.write_text(json.dumps({"theme": "presentation", "overrides": {}}))
    run(["Rscript", str(SCRIPTS / "R_mc_profile_by_context.R"), *args,
         str(work / "new_presentation.pdf"), "mean methylation", str(theme)])
    renders.append(("new (presentation theme)", work / "new_presentation.pdf"))
    return renders


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


FAMILIES = {
    "mc_profile": family_mc_profile,
    "scales": family_scales,
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
            renders = FAMILIES[fam](work, old_scripts, old_smk)
            pngs = [(t, rasterise(p)) for t, p in renders]
            side_by_side(pngs, work / "compare.png")
            old_md5 = pixels_md5(pngs[0][1])
            report = [f"family: {fam}", f"old ref: {ref}"]
            for (t, p), (_, png) in zip(renders, pngs):
                same = "identical to old" if pixels_md5(png) == old_md5 else "differs from old"
                report.append(f"{t}: {p.name}, {p.stat().st_size} bytes, pixels {same}")
            (work / "report.txt").write_text("\n".join(report) + "\n")
            print("\n".join(report), f"\n-> {work / 'compare.png'}\n")


if __name__ == "__main__":
    main()
