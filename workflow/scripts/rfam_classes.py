"""Choose which Rfam families count as structural RNA for sRNA depletion.

Rfam is not a structural-RNA database: alongside rRNA, tRNA and snRNA it
holds miRNA precursors, lncRNAs, riboswitches, CRISPR repeats and more.
Depleting every Rfam hit therefore removes reads from loci an sRNA analysis
is usually after (miRNAs above all), so the auto-derived depletion FASTA is
built from a chosen subset of Rfam classes.

Where the class comes from
--------------------------
The covariance-model file (``Rfam.cm``) carries no class annotation. Each
family's class is the ``type`` column of ``database_files/family.txt.gz``
(column 19; the same text as the ``#=GF TP`` line of ``Rfam.seed``). It is a
``;``-separated path through a small hierarchy::

    Gene; rRNA;
    Gene; tRNA;
    Gene; snRNA; snoRNA; CD-box;
    Gene; snRNA; snoRNA; HACA-box;
    Gene; snRNA; snoRNA; scaRNA;
    Gene; snRNA; splicing;
    Gene; miRNA;
    Cis-reg; riboswitch;
    Intron;

A class name selects a family when it is *any* term of that path. So
``snRNA`` takes both the spliceosomal snRNAs (``snRNA; splicing``) and,
because Rfam files snoRNA under snRNA, the C/D-box, H/ACA-box and scaRNA
families; ``snoRNA`` alone takes only the latter. ``Gene`` or ``Cis-reg`` select whole branches.

Selection modes
---------------
``structural_rna_classes`` in the options file holds exactly one of

* ``filter: [classes]`` -- deplete hits to families in these classes
  (default: rRNA, tRNA, snRNA -- all snRNAs, snoRNAs included);
* ``keep: [classes]`` -- deplete every Rfam hit *except* those in these
  classes. ``keep: []`` depletes all Rfam hits (the pre-selection behaviour).
"""

import argparse
import gzip
import sys
from collections import Counter

CONFIG_KEY = "structural_rna_classes"
MODES = ("filter", "keep")
DEFAULT_MODE = "filter"
DEFAULT_CLASSES = ("rRNA", "tRNA", "snRNA")

# Every term used in the ``type`` column of Rfam 15.1's family.txt. Config
# validation runs at startup, before family.txt has been downloaded, so it
# checks against this list; the selection step re-checks against the table
# it actually reads and warns if a requested class matches no family there.
KNOWN_CLASSES = frozenset({
    "Gene", "Cis-reg", "Intron",
    # Gene; ...
    "rRNA", "tRNA", "snRNA", "snoRNA", "CD-box", "HACA-box", "scaRNA",
    "splicing", "miRNA", "sRNA", "lncRNA", "antisense", "ribozyme",
    "CRISPR", "antitoxin",
    # Cis-reg; ...
    "riboswitch", "leader", "thermoregulator", "IRES", "frameshift_element",
})

FAMILY_TABLE_TYPE_COLUMN = 18  # 0-based index of `type` in family.txt


def parse_type(type_str):
    """Split an Rfam type string into its terms: 'Gene; rRNA;' -> ('Gene', 'rRNA')."""
    return tuple(t.strip() for t in (type_str or "").split(";") if t.strip())


def validate_structural_rna_classes(config):
    """Return a list of error messages for the structural_rna_classes option.

    An absent option is valid (the default applies).
    """
    if CONFIG_KEY not in config or config[CONFIG_KEY] is None:
        return []
    opt = config[CONFIG_KEY]
    if not isinstance(opt, dict):
        return [f"{CONFIG_KEY} must be a mapping with one key, 'filter' or "
                f"'keep', holding a list of Rfam classes; got {opt!r}."]
    errors = []
    unknown_keys = sorted(set(map(str, opt)) - set(MODES))
    if unknown_keys:
        errors.append(f"{CONFIG_KEY}: unknown key(s) {', '.join(unknown_keys)}; "
                      f"use 'filter' or 'keep'.")
    present = [m for m in MODES if m in opt]
    if len(present) == 2:
        errors.append(f"{CONFIG_KEY}: set either 'filter' or 'keep', not both.")
    elif not present and not unknown_keys:
        errors.append(f"{CONFIG_KEY}: set 'filter' or 'keep' (or remove the "
                      f"option to use the default).")
    for mode in present:
        classes = opt[mode]
        if classes is None:
            classes = []
        if isinstance(classes, str) or not isinstance(classes, (list, tuple)):
            errors.append(f"{CONFIG_KEY}.{mode} must be a list of Rfam classes; "
                          f"got {classes!r}.")
            continue
        bad = sorted({str(c) for c in classes} - KNOWN_CLASSES)
        if bad:
            errors.append(
                f"{CONFIG_KEY}.{mode}: unknown Rfam class(es) {', '.join(bad)}. "
                f"Known classes (from Rfam's family `type` field): "
                f"{', '.join(sorted(KNOWN_CLASSES, key=str.lower))}.")
        if mode == "filter" and not classes:
            errors.append(
                f"{CONFIG_KEY}.filter is empty, which would deplete nothing; "
                f"set structural_rna_depletion: false instead.")
    return errors


def selection_from_config(config):
    """(mode, sorted classes) chosen by the options file; assumes it validated."""
    opt = config.get(CONFIG_KEY) or {}
    for mode in MODES:
        if mode in opt:
            return mode, tuple(sorted({str(c) for c in (opt[mode] or [])}))
    return DEFAULT_MODE, tuple(sorted(DEFAULT_CLASSES))


def selection_label(mode, classes):
    """Stable one-line form of a selection, used as a rule param."""
    return f"{mode}:{','.join(classes) if classes else '-'}"


def parse_selection_label(label):
    mode, _, rest = label.partition(":")
    if mode not in MODES:
        raise ValueError(f"bad selection {label!r}: mode must be one of {MODES}")
    classes = () if rest in ("", "-") else tuple(rest.split(","))
    return mode, classes


def read_family_table(lines):
    """Map Rfam accession -> (family id, type terms) from family.txt lines."""
    families = {}
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) <= FAMILY_TABLE_TYPE_COLUMN or not fields[0].startswith("RF"):
            continue
        families[fields[0]] = (fields[1], parse_type(fields[FAMILY_TABLE_TYPE_COLUMN]))
    return families


def is_depleted(terms, mode, classes):
    """Whether a family with these type terms is depleted under a selection.

    ``terms`` is None for an accession missing from the family table: it has
    no known class, so ``filter`` never selects it and ``keep`` never spares it.
    """
    hit = bool(terms) and any(t in classes for t in terms)
    return hit if mode == "filter" else not hit


def parse_tblout_line(line):
    """Fields used from one cmscan --fmt 2 tblout data line, or None."""
    if not line.strip() or line.startswith("#"):
        return None
    f = line.split()
    # idx target accession query ... seq_from(9) seq_to(10) strand(11)
    return {
        "family": f[1], "acc": f[2], "chrom": f[3],
        "seq_from": int(f[9]), "seq_to": int(f[10]), "strand": f[11],
    }


def bed_record(hit):
    """BED6 line for a hit, named exactly as build_structural_rna_db names it."""
    chrom, strand = hit["chrom"], hit["strand"]
    if strand == "+":
        s, e = hit["seq_from"] - 1, hit["seq_to"]
    else:
        s, e = hit["seq_to"] - 1, hit["seq_from"]
    name = f"{hit['family']}::{chrom}:{s}-{e}({strand})"
    return chrom, s, e, f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}"


def select_hits(tblout_lines, families, mode, classes):
    """Split tblout hits into depleted BED lines and per-type tallies.

    Returns (bed_lines sorted by chrom/start, Counter of depleted types,
    Counter of spared types, set of accessions not in the family table).
    """
    classes = set(classes)
    records, depleted, spared, missing = [], Counter(), Counter(), set()
    for line in tblout_lines:
        hit = parse_tblout_line(line)
        if hit is None:
            continue
        fam = families.get(hit["acc"])
        if fam is None:
            missing.add(hit["acc"])
        terms = fam[1] if fam else None
        label = "; ".join(terms) if terms else "(not in family table)"
        if is_depleted(terms, mode, classes):
            depleted[label] += 1
            records.append(bed_record(hit))
        else:
            spared[label] += 1
    records.sort(key=lambda r: (r[0], r[1], r[2]))
    return [r[3] for r in records], depleted, spared, missing


def _open(path):
    # family.txt is not clean UTF-8 (free-text comment fields carry cp1252
    # bytes); only the ASCII accession/id/type columns matter here.
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def main(argv=None):
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--family-table", required=True, help="Rfam family.txt(.gz)")
    ap.add_argument("--tblout", required=True, help="cmscan --fmt 2 tblout")
    ap.add_argument("--selection", required=True,
                    help="'filter:rRNA,tRNA' or 'keep:miRNA' (see selection_label)")
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args(argv)

    mode, classes = parse_selection_label(args.selection)
    with _open(args.family_table) as fh:
        families = read_family_table(fh)
    if not families:
        sys.exit(f"ERROR: no Rfam families read from {args.family_table}")

    available = {t for _, terms in families.values() for t in terms}
    for c in classes:
        if c not in available:
            print(f"WARNING: class '{c}' matches no family in this Rfam "
                  f"release's family table", file=sys.stderr)

    with open(args.tblout) as fh:
        bed, depleted, spared, missing = select_hits(fh, families, mode, classes)
    with open(args.out_bed, "w") as out:
        for line in bed:
            out.write(line + "\n")

    what = ", ".join(classes) if classes else "(none)"
    print(f"Selection: {mode} {what}", file=sys.stderr)
    print(f"Depleting {sum(depleted.values())} hits; sparing "
          f"{sum(spared.values())}", file=sys.stderr)
    for title, tally in (("depleted", depleted), ("spared", spared)):
        for label, n in tally.most_common():
            print(f"  {title}\t{n}\t{label}", file=sys.stderr)
    if missing:
        print(f"WARNING: {len(missing)} accession(s) hit by cmscan are absent "
              f"from the family table (Rfam.cm and family.txt from different "
              f"releases?): {', '.join(sorted(missing))}", file=sys.stderr)


if __name__ == "__main__":
    main()
