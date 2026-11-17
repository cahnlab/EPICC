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

A class name matches a family when it is *any* term of that path. So
``snRNA`` matches both the spliceosomal snRNAs (``snRNA; splicing``) and,
because Rfam files snoRNA under snRNA, the C/D-box, H/ACA-box and scaRNA
families; ``snoRNA`` alone matches only the latter. ``Gene``, ``Cis-reg`` and
``Intron`` match whole branches.

Selection
---------
``structural_rna_classes`` in the options file has two lists:

* ``filter`` -- classes to deplete (default ``[rRNA, tRNA, snRNA]``: rRNA,
  tRNA and all snRNAs, snoRNAs included). ``filter: []`` turns structural
  RNA depletion off altogether: no Infernal scan and no depletion step, even
  for a genome with an explicit ``structural_rna_fafile``.
* ``keep`` -- classes exempt from ``filter`` (default ``[]``).

A family is depleted when one of its type terms is in ``filter`` and none is
in ``keep``: ``keep`` wins. A class listed in both is therefore never depleted
(a startup warning names it). Either key may be omitted; an omitted ``filter``
takes the default. Everything except miRNAs, for example, is
``filter: [Gene, Cis-reg, Intron]`` with ``keep: [miRNA]``.

``structural_rna_depletion``, the boolean this option replaced, is still
read: ``false`` means ``filter: []`` and ``true`` is ignored, each with a
deprecation warning.
"""

import argparse
import gzip
import sys
from collections import Counter

CONFIG_KEY = "structural_rna_classes"
LEGACY_TOGGLE_KEY = "structural_rna_depletion"
LIST_KEYS = ("filter", "keep")
DEFAULT_FILTER = ("rRNA", "tRNA", "snRNA")
DEFAULT_KEEP = ()

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


def _option(config):
    opt = config.get(CONFIG_KEY)
    return opt if isinstance(opt, dict) else {}


def validate_structural_rna_classes(config):
    """Return a list of error messages for the structural_rna_classes option.

    An absent option, or one with either list omitted, is valid (defaults
    apply).
    """
    errors = []
    legacy = config.get(LEGACY_TOGGLE_KEY)
    if legacy is not None and not isinstance(legacy, bool):
        errors.append(f"{LEGACY_TOGGLE_KEY} (deprecated) must be true or "
                      f"false; got {legacy!r}. Prefer {CONFIG_KEY}.filter: [] "
                      f"to turn depletion off.")
    if config.get(CONFIG_KEY) is None:
        return errors
    opt = config[CONFIG_KEY]
    if not isinstance(opt, dict):
        return errors + [
            f"{CONFIG_KEY} must be a mapping with 'filter' and/or 'keep', "
            f"each a list of Rfam classes; got {opt!r}."]
    unknown_keys = sorted(set(map(str, opt)) - set(LIST_KEYS))
    if unknown_keys:
        errors.append(f"{CONFIG_KEY}: unknown key(s) {', '.join(unknown_keys)}; "
                      f"use 'filter' and/or 'keep'.")
    for key in LIST_KEYS:
        if key not in opt:
            continue
        classes = opt[key]
        if isinstance(classes, str) or not isinstance(classes, (list, tuple)):
            hint = " (use [] for an empty list)" if classes is None else ""
            errors.append(f"{CONFIG_KEY}.{key} must be a list of Rfam classes"
                          f"{hint}; got {classes!r}.")
            continue
        bad = sorted({str(c) for c in classes} - KNOWN_CLASSES)
        if bad:
            errors.append(
                f"{CONFIG_KEY}.{key}: unknown Rfam class(es) {', '.join(bad)}. "
                f"Known classes (from Rfam's family `type` field): "
                f"{', '.join(sorted(KNOWN_CLASSES, key=str.lower))}.")
    return errors


def structural_rna_warnings(config):
    """Startup warnings: the deprecated toggle, and classes in both lists."""
    warnings = []
    legacy = config.get(LEGACY_TOGGLE_KEY)
    if legacy is False:
        warnings.append(
            f"{LEGACY_TOGGLE_KEY}: false is deprecated; it is treated as "
            f"{CONFIG_KEY}.filter: [] (no structural RNA depletion). Set that "
            f"instead.")
    elif legacy is True:
        warnings.append(
            f"{LEGACY_TOGGLE_KEY}: true is deprecated and ignored; depletion "
            f"is on unless {CONFIG_KEY}.filter is []. Remove the line.")
    filt, keep = _raw_lists(config)
    both = sorted(set(filt) & set(keep))
    if both:
        warnings.append(
            f"{CONFIG_KEY}: {', '.join(both)} listed in both 'filter' and "
            f"'keep'; 'keep' wins, so families of "
            f"{'this class' if len(both) == 1 else 'these classes'} are not "
            f"depleted.")
    return warnings


def _raw_lists(config):
    """(filter, keep) as configured, defaults filled, legacy toggle applied."""
    opt = _option(config)
    filt = opt.get("filter", DEFAULT_FILTER)
    keep = opt.get("keep", DEFAULT_KEEP)
    filt = filt if isinstance(filt, (list, tuple)) else DEFAULT_FILTER
    keep = keep if isinstance(keep, (list, tuple)) else DEFAULT_KEEP
    if config.get(LEGACY_TOGGLE_KEY) is False:
        filt = ()
    return tuple(map(str, filt)), tuple(map(str, keep))


def selection_from_config(config):
    """Effective (filter, keep), sorted, as the selection step applies them.

    A class in both lists is dropped from ``filter``: ``keep`` exempts every
    family it matches, so it can never deplete anything. With that, two
    configurations that deplete the same families by construction give the
    same selection (and so the same rule param).
    """
    filt, keep = _raw_lists(config)
    keep = tuple(sorted(set(keep)))
    filt = tuple(sorted(set(filt) - set(keep)))
    return filt, keep


def depletion_enabled(config):
    """Whether sRNA structural RNA depletion runs at all (filter not empty)."""
    return bool(selection_from_config(config)[0])


def selection_label(filt, keep):
    """Stable one-line form of a selection, used as a rule param."""
    return (f"filter:{','.join(filt) if filt else '-'};"
            f"keep:{','.join(keep) if keep else '-'}")


def parse_selection_label(label):
    parts = dict(p.partition(":")[::2] for p in label.split(";"))
    if set(parts) != set(LIST_KEYS):
        raise ValueError(f"bad selection {label!r}: expected "
                         f"'filter:a,b;keep:c' (use '-' for an empty list)")
    def split(v):
        return () if v in ("", "-") else tuple(v.split(","))
    return split(parts["filter"]), split(parts["keep"])


def read_family_table(lines):
    """Map Rfam accession -> (family id, type terms) from family.txt lines."""
    families = {}
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) <= FAMILY_TABLE_TYPE_COLUMN or not fields[0].startswith("RF"):
            continue
        families[fields[0]] = (fields[1], parse_type(fields[FAMILY_TABLE_TYPE_COLUMN]))
    return families


def is_depleted(terms, filt, keep=()):
    """Whether a family with these type terms is depleted: a term in
    ``filter`` and none in ``keep``.

    ``terms`` is None for an accession missing from the family table: it has
    no known class, so nothing selects it.
    """
    if not terms:
        return False
    return any(t in filt for t in terms) and not any(t in keep for t in terms)


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


def select_hits(tblout_lines, families, filt, keep=()):
    """Split tblout hits into depleted BED lines and per-type tallies.

    Returns (bed_lines sorted by chrom/start, Counter of depleted types,
    Counter of spared types, set of accessions not in the family table).
    """
    filt, keep = set(filt), set(keep)
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
        if is_depleted(terms, filt, keep):
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
                    help="'filter:rRNA,tRNA;keep:-' (see selection_label)")
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args(argv)

    filt, keep = parse_selection_label(args.selection)
    with _open(args.family_table) as fh:
        families = read_family_table(fh)
    if not families:
        sys.exit(f"ERROR: no Rfam families read from {args.family_table}")

    available = {t for _, terms in families.values() for t in terms}
    for c in filt + keep:
        if c not in available:
            print(f"WARNING: class '{c}' matches no family in this Rfam "
                  f"release's family table", file=sys.stderr)

    with open(args.tblout) as fh:
        bed, depleted, spared, missing = select_hits(fh, families, filt, keep)
    with open(args.out_bed, "w") as out:
        for line in bed:
            out.write(line + "\n")

    print(f"Selection: filter {', '.join(filt) or '(none)'}; "
          f"keep {', '.join(keep) or '(none)'}", file=sys.stderr)
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
