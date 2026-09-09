"""Centralized sample-sheet parsing and metadata utilities.

Replaces the scattered parsing logic from the Snakefile (lines 62-193)
with a single importable module. All sample-sheet column names, assay
vocabularies, environment mappings, and helper functions live here.
"""

import io
import re
from collections import OrderedDict
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Peak_type sits last: it applies to only three assays, so it stays out of the
# way of the fields every row needs. Column lookup is by header name (pandas
# reads by header), so the order is ergonomic only — but keep it in step with
# COLUMNS in tools/epicc-builder.html.
NEW_COLNAMES = [
    "Sample_ID", "Assay", "Genome", "Levels", "Replicate_ID",
    "Read_files", "Read_layout", "IP_target", "Control", "Peak_type",
    "Comments",
]

# Columns that may be absent from a sheet; filled with "" when missing. Comments
# is a free-text user annotation the pipeline never reads, so it is optional
# both ways: old sheets lack it, new sheets may leave it blank.
OPTIONAL_COLNAMES = ("Peak_type", "IP_target", "Control", "Comments")

# Assay is the *experimental method*; the peak-calling type (broad/narrow) is a
# separate analytical parameter carried in the Peak_type column, NOT baked into
# the assay name. These three pulldown assays take a Peak_type. ATAC has peaks
# but a fixed (narrow) type and no pulldown, so it carries no Peak_type.
PEAK_TYPE_ASSAYS = {"ChIP", "CUT_RUN", "CUT_TAG"}
VALID_PEAK_TYPES = {"broad", "narrow"}

# Internal combined tokens, produced from (Assay, Peak_type) by
# combine_assay_peaktype and written back over the Assay column in
# add_compat_columns. Downstream naming and the env/peaktype lookups below key
# on these, so they are part of the internal vocabulary — they are NOT accepted
# in a sample sheet's Assay column.
_COMBINED_PEAK_ASSAYS = [
    "ChIP_broad", "ChIP_narrow",
    "CUT_RUN_broad", "CUT_RUN_narrow",
    "CUT_TAG_broad", "CUT_TAG_narrow",
]

# Accepted in the Assay column. Validation (check_table) runs on this list
# before add_compat_columns folds in the peak type.
VALID_ASSAYS = [
    "ChIP", "CUT_RUN", "CUT_TAG",
    "ATAC",
    "RNAseq", "RAMPAGE",
    "sRNA",
    "WGBS", "WGBS_nd", "PBAT", "EMseq", "dmC",
]

# Assays that pull down a target via antibody and call peaks against an
# (Input/WCE/IgG) control. They share the ChIP env, peak-type machinery, and
# IP_target/Control sample-sheet semantics. Holds the separated and combined
# forms both: add_compat_columns overwrites Assay with the combined token, so
# membership tests run against either form depending on the caller.
IP_PEAK_ASSAYS = {
    "ChIP", "CUT_RUN", "CUT_TAG",
    *_COMBINED_PEAK_ASSAYS,
}

ASSAY_TO_ENV = {
    "ChIP": "ChIP", "CUT_RUN": "ChIP", "CUT_TAG": "ChIP",
    "ChIP_broad": "ChIP", "ChIP_narrow": "ChIP",
    "CUT_RUN_broad": "ChIP", "CUT_RUN_narrow": "ChIP",
    "CUT_TAG_broad": "ChIP", "CUT_TAG_narrow": "ChIP",
    "ATAC": "ATAC",
    "RNAseq": "RNA",
    "RAMPAGE": "RNA",
    "sRNA": "sRNA",
    "WGBS": "mC",
    "WGBS_nd": "mC",
    "PBAT": "mC",
    "EMseq": "mC",
    "dmC": "mC",
}

# Keyed on the internal *combined* token (produced by combine_assay_peaktype).
ASSAY_TO_PEAKTYPE = {
    "ChIP_broad": "broad",
    "ChIP_narrow": "narrow",
    "CUT_RUN_broad": "broad",
    "CUT_RUN_narrow": "narrow",
    "CUT_TAG_broad": "broad",
    "CUT_TAG_narrow": "narrow",
    "ATAC": "narrow",
}


def combine_assay_peaktype(assay, peak_type):
    """Fold a separated (Assay, Peak_type) pair into the internal combined
    assay token used by downstream naming and lookups.

    - Pulldown form ('ChIP' + 'broad')          -> 'ChIP_broad'
    - Non-peak assays / ATAC / missing peak_type -> assay unchanged

    A missing peak type on a pulldown assay is left as the bare assay so
    validation surfaces a clear error rather than silently guessing. Applying
    this to an already-combined token is a no-op, since the combined form is
    not in PEAK_TYPE_ASSAYS.
    """
    assay = (assay or "").strip()
    peak_type = (peak_type or "").strip()
    if assay in PEAK_TYPE_ASSAYS and peak_type:
        return f"{assay}_{peak_type}"
    return assay

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _row_apply(df, func):
    """Like ``df.apply(func, axis=1)`` but safe on an empty DataFrame.

    ``DataFrame.apply(..., axis=1)`` returns an empty *DataFrame* (not a
    Series) when there are no rows, which cannot be assigned back to a single
    column. An empty sheet is a legitimate input (e.g. ``snakemake --unlock``
    parses the workflow without a real sample sheet), so guard against it.
    """
    if len(df) == 0:
        return pd.Series(dtype=object)
    return df.apply(func, axis=1)


# ---------------------------------------------------------------------------
# Levels helpers
# ---------------------------------------------------------------------------

def parse_levels(levels_str):
    """Parse a Levels string into an OrderedDict of factor:level pairs.

    >>> parse_levels("genotype:WT,tissue:root")
    OrderedDict([('genotype', 'WT'), ('tissue', 'root')])
    """
    result = OrderedDict()
    for pair in levels_str.split(","):
        pair = pair.strip()
        if ":" not in pair:
            raise ValueError(
                f"Invalid Levels entry '{pair}' — expected 'factor:level' format"
            )
        factor, level = pair.split(":", 1)
        result[factor.strip()] = level.strip()
    return result


def levels_to_label(levels_str):
    """Convert Levels string to a filesystem-safe label using level values.

    >>> levels_to_label("genotype:WT,tissue:root")
    'WT_root'
    """
    parsed = parse_levels(levels_str)
    return "_".join(parsed.values())


def levels_to_factors(levels_str):
    """Return the list of factor names from a Levels string.

    >>> levels_to_factors("genotype:WT,tissue:root")
    ['genotype', 'tissue']
    """
    return list(parse_levels(levels_str).keys())


# ---------------------------------------------------------------------------
# Read-file helpers
# ---------------------------------------------------------------------------

def parse_read_files(read_files_str, read_layout):
    """Parse the Read_files field into structured information.

    Returns (components, is_sra) where:
    - components: list of individual read-file entries (after splitting on '+')
      Each component is a string (SRA ID or file path, possibly comma-separated for PE)
    - is_sra: True if the entries are SRA accessions

    For SRA with '+' merge: "SRR111+SRR222" → (["SRR111", "SRR222"], True)
    For local PE: "/a/r1.fq,/a/r2.fq" → (["/a/r1.fq,/a/r2.fq"], False)
    For local SE: "/a/reads.fq" → (["/a/reads.fq"], False)
    """
    parts = [p.strip() for p in read_files_str.split("+")]
    # Check if any part looks like an SRA/ENA/DDBJ accession (SRR/ERR/DRR)
    is_sra = bool(re.match(r"^[SED]RR\d+$", parts[0]))
    return parts, is_sra


_FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# Extensions that identify an input's type. Kept next to the FASTQ list so the
# two cannot drift apart.
_KNOWN_EXTENSIONS = (".bam", ".bed.gz", ".bedmethyl.gz") + _FASTQ_EXTENSIONS


def _strip_version_suffix(name):
    """Drop a trailing ``.<digits>`` version component from a filename.

    SRA's "Original Format" objects carry one — the submitted BAM for a run is
    served as e.g. ``Nvec_..._Nanopore.bam.1`` — which hides the real extension
    from the type checks below.

    Deliberately conservative: the suffix is removed ONLY when doing so exposes
    a known extension. So ``reads.bam.1`` -> ``reads.bam``, while a file
    genuinely named ``sample.1`` or ``chunk.12`` is left untouched and still
    reported as unrecognized rather than being silently mistyped.
    """
    base, sep, tail = name.rpartition(".")
    if not sep or not tail.isdigit():
        return name
    return base if base.endswith(_KNOWN_EXTENSIONS) else name


def _is_url(path):
    """Return True if the path looks like an HTTP(S) URL."""
    return path.startswith("http://") or path.startswith("https://")


def _is_s3_uri(path):
    """Return True if the path is an ``s3://bucket/key`` URI."""
    return path.startswith("s3://")


def _is_remote(path):
    """True for anything fetched over the network rather than read from disk.

    HTTP(S) URLs and s3:// URIs. Kept separate from ``_is_url`` so that helper
    keeps its narrower "is literally an http(s) URL" meaning.
    """
    return _is_url(path) or _is_s3_uri(path)


def s3_uri_to_https(uri):
    """Translate ``s3://bucket/key`` to a public HTTPS URL.

    EPICC supports *public* (authentication-free) S3 objects only: no signing,
    no credentials, no aws CLI. Such objects are plain HTTPS GETs, so the whole
    of S3 support is this string transform — every download path already knows
    how to curl an HTTPS URL.

    Region: we emit the region-less global endpoint rather than guessing a
    region. S3 answers it with a 307 redirect to the bucket's real region, and
    every curl in the workflow already passes ``--location``, so buckets in any
    region resolve. Requiring users to spell out a region would be both more
    typing and easy to get wrong.

    Addressing style: virtual-hosted (``https://bucket.s3.amazonaws.com/key``)
    normally, but path-style (``https://s3.amazonaws.com/bucket/key``) when the
    bucket name contains a dot — the wildcard cert for ``*.s3.amazonaws.com``
    does not match a dotted bucket label, so virtual-hosted style would fail
    TLS verification for those buckets.

    A non-s3:// value is returned unchanged, so this is safe to apply blindly.
    """
    if not _is_s3_uri(uri):
        return uri
    remainder = uri[len("s3://"):]
    bucket, sep, key = remainder.partition("/")
    if not bucket or not sep or not key:
        raise ValueError(
            f"Malformed S3 URI '{uri}' — expected 's3://bucket/key'"
        )
    if "." in bucket:
        return f"https://s3.amazonaws.com/{bucket}/{key}"
    return f"https://{bucket}.s3.amazonaws.com/{key}"


def _resolve_remote(token):
    """Map one Read_files token to the URL/path a download rule should use.

    s3:// URIs become their public HTTPS equivalent; everything else is
    unchanged.
    """
    return s3_uri_to_https(token) if _is_s3_uri(token) else token


def _resolve_remote_parts(parts):
    """Apply :func:`_resolve_remote` across '+'-joined, ','-paired components."""
    return "+".join(
        ",".join(_resolve_remote(mate.strip()) for mate in comp.split(","))
        for comp in parts
    )


def _url_basename(url):
    """Extract the filename from a URL, stripping query parameters."""
    import os
    path_part = url.split("?")[0].split("#")[0]
    return os.path.basename(path_part)


def get_seq_id_and_path(read_files_str, read_layout):
    """Bridge function: derive old-style seq_id and fastq_path from Read_files.

    Returns (seq_id, fastq_path) for backward compatibility with existing
    download rules that use these fields.

    Read_files formats handled:
    - SRA accession(s): "SRR111" or "SRR111+SRR222" → seq_id="SRR111,SRR222", path="SRA"
    - BAM file: "/path/to/file.bam" → seq_id="file", path="/path/to/file.bam"
    - BAM URL: "https://host/file.bam" → seq_id="file", path="https://host/file.bam"
    - bedMethyl: "/path/to/file.bed.gz" → seq_id="file", path="/path/to/file.bed.gz"
    - bedMethyl URL: "https://host/file.bed.gz" → seq_id="file", path="https://..."
    - Explicit FASTQ paths: "r1.fq.gz,r2.fq.gz" → seq_id="EXPLICIT", path="r1.fq.gz,r2.fq.gz"
    - FASTQ URL(s): "https://host/r.fq.gz" → seq_id="URL", path="https://host/r.fq.gz"

    Public S3 (``s3://bucket/key``) is treated as a remote URL and translated to
    its public HTTPS endpoint here — the single boundary where the path handed
    to the download rules is produced, so those rules need no S3 awareness.
    """
    parts, is_sra = parse_read_files(read_files_str, read_layout)
    if is_sra:
        # SRA: seq_id is the comma-joined accessions, path is "SRA"
        return ",".join(parts), "SRA"
    else:
        import os
        first_file = parts[0].split(",")[0]  # first mate for PE
        # Strip query params for extension checks on remote inputs (URLs and
        # s3:// URIs alike — an S3 key can carry a query string too), then drop
        # any SRA version suffix so the real extension is visible.
        check_name = _url_basename(first_file) if _is_remote(first_file) else first_file
        check_name = _strip_version_suffix(check_name)
        if check_name.endswith(".bam"):
            # BAM input (local path or URL): path IS the file/URL
            seq_id = os.path.splitext(os.path.basename(check_name))[0]
            return seq_id, _resolve_remote(first_file)
        elif check_name.endswith(".bed.gz") or check_name.endswith(".bedmethyl.gz"):
            # bedMethyl input (local path or URL): path IS the file/URL
            base = os.path.basename(check_name)
            for suffix in (".bedmethyl.gz", ".bed.gz"):
                if base.endswith(suffix):
                    seq_id = base[:-len(suffix)]
                    break
            return seq_id, _resolve_remote(first_file)
        elif any(check_name.endswith(ext) for ext in _FASTQ_EXTENSIONS):
            # Carry ALL '+'-components so the download rules can merge them
            # (each component is one SE file or a comma-separated R1,R2 pair).
            if _is_remote(first_file):
                # FASTQ URL(s): use "URL" sentinel so download rules can dispatch
                return "URL", _resolve_remote_parts(parts)
            else:
                # Explicit FASTQ path(s): pass through the full Read_files string
                return "EXPLICIT", "+".join(parts)
        else:
            raise ValueError(
                f"Unrecognized Read_files format: '{read_files_str}'. "
                f"Expected SRA accession, .bam, .bed.gz, or explicit FASTQ "
                f"path ending in .fastq.gz/.fq.gz/.fastq/.fq (as a local path, "
                f"an HTTP(S) URL, or a public s3://bucket/key URI)"
            )


# ---------------------------------------------------------------------------
# Sample sheet I/O
# ---------------------------------------------------------------------------

def read_sample_sheet(filepath):
    """Read and validate a new-format sample sheet TSV.

    Full-line comments (lines whose first non-whitespace character is '#')
    are skipped — useful for parking samples without deleting their rows.
    A '#' *within* a field (e.g. a URL fragment or a path) is preserved.

    Returns a DataFrame sorted by [Genome, Assay, Levels, Sample_ID].
    """
    # Drop full-line comments before parsing. Note: we deliberately do NOT
    # use pandas' comment="#", which truncates every line at the first '#'
    # and would corrupt any Read_files URL/path containing one.
    with open(filepath) as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", header=0, dtype=str)
    df.columns = df.columns.str.strip()

    # Fill optional columns with empty strings: Peak_type applies to only the
    # three pulldown assays and Comments to none, so a sheet may omit either
    # column outright (e.g. an all-WGBS sheet needs neither).
    for col in NEW_COLNAMES:
        if col not in df.columns:
            if col in OPTIONAL_COLNAMES:
                df[col] = ""
            else:
                raise ValueError(f"Required column '{col}' missing from sample sheet")

    # Normalize NaN to empty string for optional columns
    for col in OPTIONAL_COLNAMES:
        df[col] = df[col].fillna("")

    df = explode_genomes(df)

    df = df.sort_values(
        by=["Genome", "Assay", "Levels", "Sample_ID"]
    ).reset_index(drop=True)

    return df


def parse_genomes(genome_str):
    """Split a Genome cell into its list of reference genomes.

    ``Genome`` accepts a comma-separated list so one sample can be mapped to
    several references (e.g. ``B73,W22``) — the reads are identical, only the
    alignment target differs. A single genome is just a one-element list.
    """
    return [g.strip() for g in str(genome_str or "").split(",") if g.strip()]


def explode_genomes(df):
    """One input row per sample becomes one internal row per (sample, genome).

    Downstream code — env mapping, analysis keys, every rule — assumes exactly
    one reference genome per row. Exploding here keeps that assumption intact
    while letting the sheet express multi-genome mapping on a single line, so
    the reads are declared (and fetched, and trimmed) only once.

    ``Sample_ID`` is deliberately NOT made unique per genome: it identifies the
    library, not the alignment. Post-alignment paths are disambiguated by the
    separate ``mapped_name`` token (see add_compat_columns).
    """
    if not len(df):
        return df
    if not df["Genome"].astype(str).str.contains(",").any():
        return df  # fast path: nothing to explode
    df = df.copy()
    df["Genome"] = df["Genome"].apply(parse_genomes)
    # Duplicates must be caught here: after the explode they are indistinguishable
    # from two legitimately separate rows, and would silently produce a doubled
    # DAG for the same (sample, genome).
    for _, row in df.iterrows():
        names = row["Genome"]
        if len(set(names)) != len(names):
            raise ValueError(
                f"Sample '{row['Sample_ID']}': Genome lists the same reference "
                f"more than once ({', '.join(names)})"
            )
    return df.explode("Genome", ignore_index=True)


# ---------------------------------------------------------------------------
# Analysis key / name building
# ---------------------------------------------------------------------------

def build_analysis_key(row):
    """Build an analysis-level grouping key from a DataFrame row.

    Returns a tuple: (Assay, levels_label, IP_target, Genome)
    """
    levels_label = levels_to_label(row["Levels"])
    ip_target = row.get("IP_target", "") or ""
    return (row["Assay"], levels_label, ip_target, row["Genome"])


def build_analysis_name(row):
    """Build an analysis-level name string from a DataFrame row.

    Components are joined with '__', omitting empty parts (e.g. blank IP_target).
    Examples:
      ChIP: "ChIP_broad__WT_root__H3K9me2__ColCEN"
      RNA:  "RNAseq__WT_root__ColCEN"
    """
    key = build_analysis_key(row)
    return "__".join(part for part in key if part)


def build_control_merge_key(row):
    """Build a merge-key for control samples that ignores Assay.

    Returns a tuple: (levels_label, IP_target, Genome). Used to group
    control-sample rows (Input/WCE/IgG) into one merged-replicate group
    when they represent the same biological control material, regardless
    of which Assay the rows happen to be labeled with. A single Input
    that serves as control for both broad and narrow ChIP IPs merges
    correctly as long as Levels, IP_target, and Genome match — even if
    rep1 is labeled ChIP_broad and rep2 ChIP_narrow.
    """
    levels_label = levels_to_label(row["Levels"])
    ip_target = row.get("IP_target", "") or ""
    return (levels_label, ip_target, row["Genome"])


def build_analysis_to_replicates(df):
    """Build a dict mapping analysis_key → list of Replicate_IDs.

    Includes every row that can be analysed in its own right, i.e. every row
    ``is_peak_call_target`` accepts. A sample that serves as another row's
    control is NOT excluded: as long as it declares a ``Control`` of its own it
    is a legitimate analysis target too (dual-role samples — e.g. an H3 ChIP
    used both as its own target and as the control for H3K9me2).
    """
    non_control = peak_callable_rows(df)
    result = {}
    for _, row in non_control.iterrows():
        key = build_analysis_key(row)
        result.setdefault(key, []).append(row["Replicate_ID"])
    return result


def get_replicate_sample_ids(analysis_name, df):
    """Get the genome-qualified replicate names for an analysis group.

    Returns ``mapped_name`` values ('{Sample_ID}__{Genome}'), not bare
    Sample_IDs: every caller uses these to build post-alignment paths, where the
    reference is part of the file's identity. The *lookup key* is still a bare
    analysis name or control Sample_ID.


    The analysis_name can be:
    - An analysis name like "ChIP_broad__WT_cell__H3K9me2__Spombe"
    - A control Sample_ID like "WT_cell_Input_rep1" (for merging controls)

    Returns the list of Sample_IDs that belong to that analysis group.
    """
    controls = identify_control_samples(df)
    # Analysable rows (dual-role samples included — see build_analysis_to_replicates).
    # A dual-role sample is reachable BOTH ways: by its analysis name here, and
    # by its Sample_ID through the control-merge fallback below. The two lookups
    # use different key formats, so they cannot collide.
    non_control = peak_callable_rows(df)

    # Try matching non-control analysis names first
    result = []
    for _, row in non_control.iterrows():
        if build_analysis_name(row) == analysis_name:
            result.append(row.get("mapped_name", row["Sample_ID"]))
    if result:
        return result

    # If no match, check whether the key names a control sample and collect
    # every control in the same merge group. build_control_merge_key is
    # Assay-agnostic, so one biological Input/IgG/WCE merges across reps whose
    # Levels/IP_target/Genome match regardless of each row's Assay label.
    #
    # match_sample_rows, not a bare ==: the key arrives genome-qualified here
    # (merging is post-alignment) while identify_control_samples yields bare
    # Sample_IDs. It also pins one (Sample_ID, Genome) row, and the merge key
    # includes Genome, so multi-genome sheets never mix references.
    control_df = df[df["Sample_ID"].isin(controls)]
    match = match_sample_rows(control_df, analysis_name)
    if not match.empty:
        ctrl_key = build_control_merge_key(match.iloc[0])
        for _, row in control_df.iterrows():
            if build_control_merge_key(row) == ctrl_key:
                result.append(row.get("mapped_name", row["Sample_ID"]))

    return result


# ---------------------------------------------------------------------------
# Control handling
# ---------------------------------------------------------------------------

def identify_control_samples(df):
    """Return the set of Sample_IDs that are referenced as controls.

    Note: this only finds controls that some *other* row points at. Use
    ``is_peak_call_target`` to decide whether a row can be peak-called — a
    control that no IP references ("orphan control") is absent from this set.
    """
    controls = set()
    for val in df["Control"]:
        if pd.notna(val) and str(val).strip():
            controls.add(str(val).strip())
    return controls


def match_sample_rows(df, name):
    """Rows for a sample named either bare or genome-qualified.

    Pre-alignment code passes a bare ``Sample_ID``; post-alignment code passes
    the ``mapped_name`` ('{Sample_ID}__{Genome}'). Accepting both here means the
    lookup helpers work on either side of the alignment boundary without every
    caller having to know which token it holds.
    """
    match = df.loc[df["Sample_ID"] == name]
    if match.empty and "mapped_name" in df.columns:
        match = df.loc[df["mapped_name"] == name]
    return match


def is_peak_call_target(row):
    """True if a row may be enumerated as a peak-calling target.

    Peak calling for the pulldown assays (``IP_PEAK_ASSAYS``) requires a
    control, so a pulldown row is a valid target only when it declares a
    non-empty ``Control``. This is a *structural* test, unlike
    ``identify_control_samples`` which only recognizes controls that some
    other row references by Sample_ID. Two cases it catches that the
    reference-based test misses:

    - **Orphan controls**: an Input/WCE/IgG row that no IP points at. It has
      no ``Control`` of its own, so ``assign_chip_input`` cannot resolve one
      and peak calling raises "No control found".
    - **IPs missing a control**: equally un-peak-callable, and they would
      fail in exactly the same way.

    Rows of non-pulldown assays (ATAC, RNAseq, sRNA, mC, ...) return True:
    the ``Control`` column does not apply to them, and their own target
    enumeration handles them (ATAC calls peaks without a control).
    """
    if row.get("Assay") not in IP_PEAK_ASSAYS:
        return True
    control = row.get("Control")
    if control is None or pd.isna(control):
        return False
    return bool(str(control).strip())


def peak_callable_rows(df):
    """Return the subset of ``df`` that ``is_peak_call_target`` accepts.

    Single place for the row-wise mask so every caller agrees on which rows are
    analysable. Note ``.apply(axis=1)`` returns a DataFrame rather than a Series
    on an empty frame, so guard the empty case.
    """
    if not len(df):
        return df.copy()
    return df[_row_apply(df, is_peak_call_target).astype(bool)].copy()


def get_control_sample_id(sample_id, df):
    """Get the control for a given sample.

    The returned name matches the form of the input: asked with a
    genome-qualified name, the control comes back qualified with the SAME
    genome, because the caller is building a post-alignment path (the control's
    BAM) and a control aligned to a different reference would be wrong.
    """
    match = match_sample_rows(df, sample_id)
    if match.empty:
        return None
    row = match.iloc[0]
    control = str(row["Control"]).strip()
    if not control or control == "nan":
        return None
    # Qualify the control iff the caller used a qualified name.
    if "mapped_name" in df.columns and sample_id == row.get("mapped_name"):
        return f"{control}__{row['ref_genome']}"
    return control


# ---------------------------------------------------------------------------
# Analysis samples (non-control)
# ---------------------------------------------------------------------------

def get_analysis_samples(df):
    """Return a DataFrame of analysis-level samples (controls filtered out).

    Deduplicated by analysis key (one row per unique analysis group).
    When an analysis group has mixed PE/SE replicates, the analysis-level
    Read_layout is forced to SE and a warning is printed (PE reads are
    treated as SE for merged peak calling to avoid MACS2 errors).
    """
    import sys

    # Keep every row that can be peak-called: a pulldown row must declare a
    # Control (this drops orphan controls and IPs with no Control, which would
    # otherwise earn selected_peaks/IDR targets they can never satisfy). Being
    # referenced as someone else's control is NOT disqualifying — dual-role
    # samples get the full analysis treatment. See is_peak_call_target.
    non_control = peak_callable_rows(df)

    # Deduplicate by analysis key
    non_control["_analysis_key"] = _row_apply(
        non_control, lambda row: build_analysis_key(row)
    )

    # Detect mixed PE/SE within analysis groups before deduplication
    mixed_groups = []
    for key, group in non_control.groupby("_analysis_key"):
        layouts = group["Read_layout"].unique()
        if len(layouts) > 1:
            name = build_analysis_name(group.iloc[0])
            sids = group["Sample_ID"].tolist()
            mixed_groups.append((key, name, sids))

    deduped = non_control.drop_duplicates(subset=["_analysis_key"]).copy()

    # Force SE for mixed groups (PE reads treated as SE for merged analysis)
    if mixed_groups:
        for key, name, sids in mixed_groups:
            mask = deduped["_analysis_key"] == key
            deduped.loc[mask, "Read_layout"] = "SE"
            sid_list = ", ".join(sids)
            print(
                f"WARNING: Analysis group '{name}' has mixed PE/SE "
                f"replicates ({sid_list}). Merged peak calling will treat "
                f"all reads as SE. Per-replicate peak calls use each "
                f"sample's own layout.",
                file=sys.stderr,
            )

    deduped = deduped.drop(columns=["_analysis_key"])
    return deduped


# ---------------------------------------------------------------------------
# Field lookup
# ---------------------------------------------------------------------------

def get_sample_field(sample_id, df, field):
    """Look up a single field value by bare or genome-qualified sample name."""
    match = match_sample_rows(df, sample_id)
    if match.empty:
        return None
    return match[field].iloc[0]


# ---------------------------------------------------------------------------
# Environment and peak type
# ---------------------------------------------------------------------------

def get_env(assay):
    """Map an Assay value to its pipeline environment folder name."""
    env = ASSAY_TO_ENV.get(assay)
    if env is None:
        raise ValueError(f"Unknown assay '{assay}' — not in ASSAY_TO_ENV")
    return env


def get_peaktype(assay, config_override=None):
    """Get peak type (broad/narrow) for a given Assay.

    If config_override is provided (a dict mapping assay→peaktype),
    it takes precedence.
    """
    if config_override and assay in config_override:
        return config_override[assay]
    pt = ASSAY_TO_PEAKTYPE.get(assay)
    if pt is None:
        raise ValueError(
            f"No peak type defined for assay '{assay}'. "
            "Peak types are only defined for ChIP_broad, ChIP_narrow, "
            "CUT_RUN_broad, CUT_RUN_narrow, CUT_TAG_broad, "
            "CUT_TAG_narrow, and ATAC."
        )
    return pt


# ---------------------------------------------------------------------------
# Backward-compatibility bridge columns
# ---------------------------------------------------------------------------

def add_compat_columns(df):
    """Add backward-compatible columns derived from new-format fields.

    These columns bridge the gap during rule migration, allowing rules
    that still reference old column names to work with new-format data.

    Adds: data_type, ref_genome, replicate, paired, line, tissue,
          sample_type, extra_info, levels_label, env, sample_name
    """
    df = df.copy()

    # Fold (Assay, Peak_type) into the internal *combined* assay token (e.g.
    # 'ChIP' + 'broad' -> 'ChIP_broad') that every downstream rule, analysis
    # name, and env/peaktype lookup keys on. check_table runs on the separated
    # form BEFORE this, so bad Peak_type values are reported there rather than
    # silently combined here.
    if "Peak_type" not in df.columns:
        df["Peak_type"] = ""
    df["Assay"] = _row_apply(
        df, lambda r: combine_assay_peaktype(r["Assay"], r.get("Peak_type", ""))
    ) if len(df) else df["Assay"]

    # Direct mappings
    df["data_type"] = df["Assay"]
    df["ref_genome"] = df["Genome"]
    df["replicate"] = df["Replicate_ID"]
    df["paired"] = df["Read_layout"]

    # Levels → line (1st factor value) and tissue (2nd factor value)
    def _extract_level(levels_str, index):
        parsed = parse_levels(levels_str)
        vals = list(parsed.values())
        if index < len(vals):
            return vals[index]
        return ""

    df["line"] = df["Levels"].apply(lambda x: _extract_level(x, 0))
    df["tissue"] = df["Levels"].apply(lambda x: _extract_level(x, 1))

    # sample_type: IP_target for IP-with-peaks assays, else Assay
    def _derive_sample_type(row):
        if row["Assay"] in IP_PEAK_ASSAYS:
            return row["IP_target"] if row["IP_target"] else row["Assay"]
        return row["Assay"]

    df["sample_type"] = _row_apply(df, _derive_sample_type)

    # extra_info: IP_target or "N/A"
    df["extra_info"] = df["IP_target"].apply(
        lambda x: x if x else "N/A"
    )

    # levels_label
    df["levels_label"] = df["Levels"].apply(levels_to_label)

    # env
    df["env"] = df["Assay"].map(ASSAY_TO_ENV)

    # sample_name = Sample_ID (the key bridge)
    df["sample_name"] = df["Sample_ID"]

    # Post-alignment identity: pre-alignment paths keep the bare sample_name
    # (one download+trim serves every reference), while alignment onward uses
    # mapped_name so two genomes cannot collide on one filename (#39). '__' is
    # the reserved delimiter and both halves are validated to exclude it, so
    # rsplit('__', 1) recovers the pair unambiguously.
    df["mapped_name"] = df["Sample_ID"].astype(str) + "__" + df["ref_genome"].astype(str)

    # Derive seq_id and fastq_path from Read_files
    def _derive_seq_id(row):
        seq_id, _ = get_seq_id_and_path(row["Read_files"], row["Read_layout"])
        return seq_id

    def _derive_fastq_path(row):
        _, fq_path = get_seq_id_and_path(row["Read_files"], row["Read_layout"])
        return fq_path

    df["seq_id"] = _row_apply(df, _derive_seq_id)
    df["fastq_path"] = _row_apply(df, _derive_fastq_path)

    return df
