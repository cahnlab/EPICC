#!/usr/bin/env python3
"""Build the `all` env's matrix by joining the `most` and mC matrices.

The `all` env holds every data type in one matrix. It used to be computed
directly from every bigwig with the non-mC computeMatrix options
(`--missingDataAsZero --skipZeros`, bin size `bs`), so its mC columns read an
uncovered cytosine as 0 % methylation instead of missing, and used the non-mC
bin size. The `most` matrix (non-mC, built with those options) and the mC
matrix (built with `base_mc`/`bs_mc`, which keep NaN) are each right for their
own data, so `all` is their join rather than a third computeMatrix run:

- columns: the `most` columns, then the mC columns.
- bins: one bin size for every column, the greatest common divisor of the two
  (`bs` 50 and `bs_mc` 100 by default give 50), each bin repeated to fill it.
  That is exact -- every value is the one computeMatrix gave, at its own
  resolution -- and it keeps one bin layout across the matrix. deeptools'
  `plotProfile --perGroup` and the R profiles (matrix_io.R) draw every line on
  one bin-index axis with the first sample's ticks, so mixed bin sizes would
  squeeze the mC lines into part of the panel. computeMatrix requires every
  segment (upstream, body, downstream, unscaled ends) to be a multiple of the
  bin size, so it is a multiple of the common divisor too.
- rows: per region group, joined on (chrom, start, end, name, strand), a
  duplicated key matched in order. computeMatrixOperations cbind is not used:
  it joins on the region name alone, collapsing duplicate names, and keeps only
  the first matrix's rows.
- row order: the mC matrix's, which has no --skipZeros and so carries every
  region in target-file order; a `most` row it lacks is appended to its group.
- fill: a region missing from `most` was dropped by --skipZeros, i.e. every
  non-mC value was zero (missing data is zero there), so it is filled with 0;
  a region missing from the mC matrix is filled with NaN.
- a row with no data in any column -- every non-mC value zero and every mC
  value NaN -- is dropped, the rule --skipZeros applies within `most`.

Either input may be omitted (an analysis with no mC, or nothing but mC); the
matrix is then rewritten from the one given, with the same drop rule.

Usage: matrix_join_all.py -o OUT.gz [--most MOST.gz] [--mc MC.gz]
"""

import argparse
import gzip
import json
import math
from collections import defaultdict, deque

# Header entries deeptools keeps one per sample column (heatmapper.special_params).
SPECIAL_PARAMS = ("unscaled 5 prime", "unscaled 3 prime", "body", "downstream",
                  "upstream", "ref point", "bin size")
# The bin layout in base pairs; must agree between the two matrices.
LAYOUT_PARAMS = ("upstream", "downstream", "body", "unscaled 5 prime",
                 "unscaled 3 prime", "ref point")


def read_matrix(path):
    """(header, {group: [(key, bed_fields, values)]}) of a computeMatrix file.

    bed_fields are the first six columns as written (chrom, start, end, name,
    score, strand); values are the remaining columns, kept as strings so the
    output reproduces them exactly.
    """
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline().lstrip("@"))
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    bounds = header["group_boundaries"]
    groups = {}
    for i, label in enumerate(header["group_labels"]):
        groups[label] = [
            ((r[0], r[1], r[2], r[3], r[5]), r[:6], r[6:])
            for r in rows[bounds[i]:bounds[i + 1]]
        ]
    return header, groups


def _per_sample(header, key):
    n = len(header["sample_labels"])
    v = header.get(key)
    return list(v) if isinstance(v, list) else [v] * n


def _resampler(header, target):
    """(new sample boundaries, fn(values) -> values at bin size `target`)."""
    sizes = _per_sample(header, "bin size")
    bounds = header["sample_boundaries"]
    reps, new_bounds = [], [0]
    for i, size in enumerate(sizes):
        if size % target:
            raise ValueError(f"bin size {size} is not a multiple of {target}")
        reps.append(size // target)
        new_bounds.append(new_bounds[-1] + (bounds[i + 1] - bounds[i]) * reps[-1])
    if all(r == 1 for r in reps):
        return new_bounds, lambda values: values

    def resample(values):
        out = []
        for i, r in enumerate(reps):
            for v in values[bounds[i]:bounds[i + 1]]:
                out.extend([v] * r)
        return out
    return new_bounds, resample


def _is_zero_or_nan(value):
    v = float(value)
    return math.isnan(v) or v == 0.0


def _is_nan(value):
    return math.isnan(float(value))


def join(most, mc):
    """Join (header, groups) pairs; either may be None. Returns (header, rows)."""
    parts = [p for p in (most, mc) if p is not None]
    if not parts:
        raise ValueError("matrix_join_all: no input matrix")
    headers = [h for h, _ in parts]

    for key in LAYOUT_PARAMS:
        values = {json.dumps(v) for h in headers for v in _per_sample(h, key)}
        if len(values) > 1:
            raise ValueError(f"matrix_join_all: the matrices disagree on '{key}' "
                             f"({', '.join(sorted(values))}); only the bin size may differ")
    target = math.gcd(*(int(s) for h in headers for s in _per_sample(h, "bin size")))

    # Header: the first input's, sample entries extended with the second's,
    # every column at the common bin size.
    header = json.loads(json.dumps(headers[0]))
    header["sample_labels"] = [s for h in headers for s in h["sample_labels"]]
    for key in SPECIAL_PARAMS:
        header[key] = [v for h in headers for v in _per_sample(h, key)]
    header["bin size"] = [target] * len(header["sample_labels"])
    if len(parts) == 2:
        # Each half followed its own rule (the per-column fill rules above);
        # neither scalar describes the joined matrix.
        for key in ("missing data as zero", "skip zeros"):
            if key in header:
                header[key] = None

    resamplers = {}
    boundaries = [0]
    for name, part in (("most", most), ("mc", mc)):
        if part is None:
            resamplers[name] = (0, lambda values: values)
            continue
        bounds, fn = _resampler(part[0], target)
        resamplers[name] = (bounds[-1], fn)
        offset = boundaries[-1]
        boundaries.extend(b + offset for b in bounds[1:])
    header["sample_boundaries"] = boundaries
    most_cols, most_fn = resamplers["most"]
    mc_cols, mc_fn = resamplers["mc"]

    most_groups = most[1] if most else {}
    mc_groups = mc[1] if mc else {}

    # Group order: mC's (its rows are the complete set), then any only in most.
    labels = list(mc_groups) + [g for g in most_groups if g not in mc_groups]

    out_rows = []
    group_bounds = [0]
    kept_labels = []
    for label in labels:
        queues = defaultdict(deque)
        for key, bed, values in most_groups.get(label, []):
            queues[key].append((bed, most_fn(values)))
        joined = []
        for key, bed, values in mc_groups.get(label, []):
            if queues[key]:
                _, most_values = queues[key].popleft()
            else:
                most_values = ["0.000000"] * most_cols
            joined.append((bed, most_values, mc_fn(values)))
        # most rows the mC matrix lacks (computeMatrix skipped the region for
        # the mC bigwigs only): keep them, mC missing.
        for queue in queues.values():
            for bed, most_values in queue:
                joined.append((bed, most_values, ["nan"] * mc_cols))
        group_rows = [
            bed + most_values + mc_values
            for bed, most_values, mc_values in joined
            if not (all(_is_zero_or_nan(v) for v in most_values)
                    and all(_is_nan(v) for v in mc_values))
        ]
        if not group_rows:
            continue
        kept_labels.append(label)
        out_rows.extend(group_rows)
        group_bounds.append(len(out_rows))

    header["group_labels"] = kept_labels
    header["group_boundaries"] = group_bounds
    return header, out_rows


def write_matrix(path, header, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header, separators=(",", ":")) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Join the most and mC matrices into the all matrix.")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--most", help="the most env's final matrix (non-mC)")
    ap.add_argument("--mc", help="the mC env's final matrix")
    args = ap.parse_args(argv)
    most = read_matrix(args.most) if args.most else None
    mc = read_matrix(args.mc) if args.mc else None
    header, rows = join(most, mc)
    if not rows:
        raise SystemExit("matrix_join_all: no region has data in any column")
    write_matrix(args.output, header, rows)


if __name__ == "__main__":
    main()
