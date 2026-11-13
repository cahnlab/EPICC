"""Name the prepared TE annotation, and the outputs built on it, by its source.

A genome's options-file TE annotation (``genomes.<G>.te_file``) is prepared
as ``<G>__TE_file.bed`` and analysed under the target name ``all_TEs``. An
alternative annotation for one run (``genomes.<G>.alt_te_file``, which
``epicc --te-file`` sets) is prepared beside it as ``<G>__TE_file_<tag>.bed``
and analysed as ``all_TEs_<tag>``, so its heatmaps, sRNA counts and browser
tracks sit next to the default ones instead of replacing them. Swapping the
source of the one shared file would make Snakemake rebuild it, and everything
downstream of it, every time a run switched between the two.

The tag comes from the file name (``EDTA_v2.gff3.gz`` -> ``EDTA_v2``). The
alternative rule writes the source path to ``<G>__TE_file_<tag>.source``; a
later run whose different file would get the same tag gets the tag plus a
short hash of its path instead, so the two never share a file.
"""

import hashlib
import os
import re

DEFAULT_TE_TARGET = "all_TEs"
ALT_TE_KEY = "alt_te_file"

_TAG_RE = re.compile(r"[A-Za-z0-9_-]+")


def te_tag_from_path(path):
    """'/data/EDTA_v2.gff3.gz' -> 'EDTA_v2'; unsafe characters become '-'."""
    name = os.path.basename(str(path).split("?", 1)[0].rstrip("/"))
    while True:
        stem, ext = os.path.splitext(name)
        if not ext or stem == "":
            break
        name = stem
    tag = re.sub(r"[^A-Za-z0-9_-]+", "-", name)
    # '__' separates fields in output names; never let the tag contain it.
    tag = re.sub(r"_+", "_", tag).strip("_-")
    return tag or "alt"


def te_source_record(genome_dir, genome, tag):
    return os.path.join(genome_dir, genome, f"{genome}__TE_file_{tag}.source")


def resolve_te_tag(path, genome_dir, genome):
    """Tag for ``path``, made unique against an earlier alternative TE file.

    The record of an earlier run with the same tag but a different source
    means the name is taken: append 8 hex digits of the path's SHA-1.
    """
    tag = te_tag_from_path(path)
    record = te_source_record(genome_dir, genome, tag)
    try:
        with open(record) as fh:
            previous = fh.read().strip()
    except OSError:
        return tag
    if previous and previous != str(path):
        tag = f"{tag}-{hashlib.sha1(str(path).encode()).hexdigest()[:8]}"
    return tag


def te_suffixes(genomes_cfg, genome_dir):
    """``{genome: suffix}``: '' for the options-file TE file, '_<tag>' otherwise.

    An ``alt_te_file`` identical to ``te_file`` is the default annotation and
    gets no tag, so passing the options-file path on the CLI redoes nothing.
    """
    out = {}
    for genome, gcfg in (genomes_cfg or {}).items():
        alt = (gcfg or {}).get(ALT_TE_KEY)
        if alt and str(alt) != str((gcfg or {}).get("te_file") or ""):
            out[genome] = "_" + resolve_te_tag(alt, genome_dir, genome)
        else:
            out[genome] = ""
    return out


def te_source(gcfg, suffix):
    """The annotation the run uses: the alternative file when one is tagged."""
    gcfg = gcfg or {}
    return gcfg.get(ALT_TE_KEY) if suffix else gcfg.get("te_file")


def te_bed_path(genome_dir, genome, suffix):
    return os.path.join(genome_dir, genome, f"{genome}__TE_file{suffix}.bed")


def te_target_name(suffix):
    return DEFAULT_TE_TARGET + suffix


def te_suffix_from_target(target_name):
    """'all_TEs' -> '', 'all_TEs_EDTA' -> '_EDTA'; None if not a TE target."""
    if target_name == DEFAULT_TE_TARGET:
        return ""
    prefix = DEFAULT_TE_TARGET + "_"
    if target_name.startswith(prefix) and _TAG_RE.fullmatch(target_name[len(prefix):]):
        return target_name[len(DEFAULT_TE_TARGET):]
    return None
