#!/usr/bin/env python3
"""Build the `all` env's matrix by joining the `most` and mC matrices.

The `all` env holds every data type in one matrix. It used to be computed
directly from every bigwig with the non-mC computeMatrix options
(`--missingDataAsZero --skipZeros`, bin size `bs`), so its mC columns read an
uncovered cytosine as 0 % methylation instead of missing, and used the non-mC
bin size. The `most` matrix (non-mC, built with those options) and the mC
matrix (built with `base_mc`/`bs_mc`, which keep NaN) are each right for their
own data, so `all` is their join rather than a third computeMatrix run. It is
the one implementation of that join: the combined heatmap (phase 3 of
dev/docs/r-plotting-plan.md) reads this matrix rather than joining again.

- columns: the `most` columns, then the mC columns.
- bins: one bin size for every column, the greatest common divisor of the two
  (`bs` 50 and `bs_mc` 100 by default give 50), each bin repeated to fill it.
  That is exact -- every value is the text computeMatrix wrote, at its own
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
- row order: target-file order. The mC matrix has no --skipZeros, so its order
  is the target file's; a `most` row it lacks goes right after the row that
  precedes it in `most`, so it keeps its place under `--sortRegions keep`.
- fill: a region missing from `most` was dropped by --skipZeros, i.e. every
  non-mC value was zero (missing data is zero there), so it is filled with 0;
  a region missing from the mC matrix is filled with NaN.
- dropping: with both inputs, a row with no data in any column -- every non-mC
  value zero and every mC value NaN -- is dropped, the rule --skipZeros applies
  within `most`. With one input (an analysis with no mC, or nothing but mC)
  every row is kept, so `all` shows the same regions as that env.

Memory: values are never held for the whole matrix. A first pass reads only
the six BED columns of both files to plan the output order; a second pass
streams both, holding a row's values as text until it is written. As long as
`most` lists its regions in the mC matrix's order (computeMatrix keeps the
target file's order in both), only one `most` row is held at a time; any other
order still works, holding the out-of-order rows until they are needed.

Usage: matrix_join_all.py -o OUT.gz [--most MOST.gz] [--mc MC.gz]
"""

import argparse
import gzip
import json
import math
import os
from collections import defaultdict, deque

# Header entries deeptools keeps one per sample column (heatmapper.special_params).
SPECIAL_PARAMS = ("unscaled 5 prime", "unscaled 3 prime", "body", "downstream",
                  "upstream", "ref point", "bin size")
# The bin layout in base pairs; must agree between the two matrices.
LAYOUT_PARAMS = ("upstream", "downstream", "body", "unscaled 5 prime",
                 "unscaled 3 prime", "ref point")
# Output compression: computeMatrix's gzip level 9 costs several times the time
# for a few percent of size on matrices this large.
COMPRESSLEVEL = 6



def read_keys(path):
    """(header, [[key, ...] per group]) -- the BED columns only."""
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline().lstrip("@"))
        keys = []
        for line in fh:
            f = line.split("\t", 6)
            keys.append((f[0], f[1], f[2], f[3], f[5]))
    b = header["group_boundaries"]
    return header, [keys[b[i]:b[i + 1]] for i in range(len(header["group_labels"]))]


def iter_rows(path):
    """(bed6 text, values text) per row, the values as one unsplit string."""
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t", 6)
            yield "\t".join(f[:6]), (f[6] if len(f) > 6 else "")


def _per_sample(header, key):
    n = len(header["sample_labels"])
    v = header.get(key)
    return list(v) if isinstance(v, list) else [v] * n


def _repeat_counts(header, target):
    """Per-column repeat count to bring a matrix to bin size `target`."""
    sizes = _per_sample(header, "bin size")
    bounds = header["sample_boundaries"]
    counts = []
    for i, size in enumerate(sizes):
        if int(size) % target:
            raise ValueError(f"bin size {size} is not a multiple of {target}")
        counts.extend([int(size) // target] * (bounds[i + 1] - bounds[i]))
    return counts


def _resampler(counts):
    """fn(values text) -> values text with column j repeated counts[j] times."""
    if all(c == 1 for c in counts):
        return lambda text: text
    if len(set(counts)) == 1:
        r = counts[0]

        def uniform(text):
            fields = text.split("\t")
            out = [None] * (len(fields) * r)
            for k in range(r):
                out[k::r] = fields
            return "\t".join(out)
        return uniform

    def mixed(text):
        return "\t".join(v for v, c in zip(text.split("\t"), counts) for _ in range(c))
    return mixed


def _no_signal(text):
    """Every value zero or NaN (short-circuits on the first that is not)."""
    for v in text.split("\t"):
        if v != "nan" and float(v) != 0.0:
            return False
    return True


def _all_nan(text):
    for v in text.split("\t"):
        if v != "nan" and not math.isnan(float(v)):
            return False
    return True


def joined_header(headers, target):
    for key in LAYOUT_PARAMS:
        values = {json.dumps(v) for h in headers for v in _per_sample(h, key)}
        if len(values) > 1:
            raise ValueError(f"matrix_join_all: the matrices disagree on '{key}' "
                             f"({', '.join(sorted(values))}); only the bin size may differ")
    header = json.loads(json.dumps(headers[0]))
    header["sample_labels"] = [s for h in headers for s in h["sample_labels"]]
    for key in SPECIAL_PARAMS:
        header[key] = [v for h in headers for v in _per_sample(h, key)]
    header["bin size"] = [target] * len(header["sample_labels"])
    boundaries = [0]
    for h in headers:
        b = h["sample_boundaries"]
        for i, size in enumerate(_per_sample(h, "bin size")):
            boundaries.append(boundaries[-1] + (b[i + 1] - b[i]) * (int(size) // target))
    header["sample_boundaries"] = boundaries
    if len(headers) == 2:
        # Each half already carries its own fill in the values (0 from most's
        # --missingDataAsZero, NaN from mC), so the joined matrix asks deeptools
        # for neither. computeMatrix always writes booleans here.
        header["missing data as zero"] = False
        header["skip zeros"] = False
    return header


def plan(most_groups, mc_groups, most_labels, mc_labels):
    """Output order: [(label, [(most_idx|None, mc_idx|None), ...])].

    Indices are global row numbers in each file. Group order is mC's, then any
    group only `most` has.
    """
    most_by_label = dict(zip(most_labels, most_groups))
    mc_by_label = dict(zip(mc_labels, mc_groups))
    most_start, mc_start = {}, {}
    n = 0
    for label, keys in zip(most_labels, most_groups):
        most_start[label] = n
        n += len(keys)
    n = 0
    for label, keys in zip(mc_labels, mc_groups):
        mc_start[label] = n
        n += len(keys)

    out = []
    for label in list(mc_labels) + [g for g in most_labels if g not in mc_by_label]:
        mc_keys = mc_by_label.get(label, [])
        most_keys = most_by_label.get(label, [])
        slots = defaultdict(deque)
        for i, key in enumerate(mc_keys):
            slots[key].append(i)
        match = [None] * len(mc_keys)     # mc local idx -> most local idx
        after = defaultdict(list)         # mc local idx (-1 = group start) -> most-only rows
        anchor = -1
        for j, key in enumerate(most_keys):
            if slots[key]:
                i = slots[key].popleft()
                match[i] = j
                anchor = i
            else:
                after[anchor].append(j)
        mo = most_start.get(label, 0)
        mc = mc_start.get(label, 0)
        rows = [(mo + j, None) for j in after[-1]]
        for i in range(len(mc_keys)):
            rows.append((None if match[i] is None else mo + match[i], mc + i))
            rows.extend((mo + j, None) for j in after[i])
        out.append((label, rows))
    return out


class _Reader:
    """Rows of one file by global index, read forward; skipped rows are held."""

    def __init__(self, path):
        self.it = iter_rows(path) if path else iter(())
        self.pos = 0
        self.held = {}

    def get(self, idx):
        if idx in self.held:
            return self.held.pop(idx)
        while self.pos < idx:
            self.held[self.pos] = next(self.it)
            self.pos += 1
        self.pos += 1
        return next(self.it)


def join_files(most_path, mc_path, out_path):
    """Write the joined matrix; returns (rows written, rows dropped)."""
    if not (most_path or mc_path):
        raise ValueError("matrix_join_all: no input matrix")
    most_h, most_groups = read_keys(most_path) if most_path else (None, [])
    mc_h, mc_groups = read_keys(mc_path) if mc_path else (None, [])
    headers = [h for h in (most_h, mc_h) if h is not None]
    target = math.gcd(*(int(s) for h in headers for s in _per_sample(h, "bin size")))
    header = joined_header(headers, target)

    most_counts = _repeat_counts(most_h, target) if most_h else []
    mc_counts = _repeat_counts(mc_h, target) if mc_h else []
    most_fn, mc_fn = _resampler(most_counts), _resampler(mc_counts)
    most_fill = "\t".join(["0.000000"] * len(most_counts))
    mc_fill = "\t".join(["nan"] * len(mc_counts))
    both = most_h is not None and mc_h is not None

    order = plan(most_groups, mc_groups,
                 most_h["group_labels"] if most_h else [],
                 mc_h["group_labels"] if mc_h else [])
    most_r, mc_r = _Reader(most_path), _Reader(mc_path)

    # The header's group boundaries depend on which rows are dropped, so the
    # rows go to a gzip member of their own first; the output is the header
    # member followed by it (gzip, deeptools and `gzip -dc` read concatenated
    # members as one stream).
    body_path = f"{out_path}.body.tmp"
    labels, bounds, dropped = [], [0], 0
    try:
        with gzip.open(body_path, "wt", compresslevel=COMPRESSLEVEL) as body:
            for label, rows in order:
                n_group = 0
                for most_idx, mc_idx in rows:
                    most_bed, most_vals = most_r.get(most_idx) if most_idx is not None else (None, None)
                    mc_bed, mc_vals = mc_r.get(mc_idx) if mc_idx is not None else (None, None)
                    if both and ((most_vals is None or _no_signal(most_vals))
                                 and (mc_vals is None or _all_nan(mc_vals))):
                        dropped += 1
                        continue
                    parts = [mc_bed if mc_bed is not None else most_bed]
                    if most_counts:
                        parts.append(most_fn(most_vals) if most_vals is not None else most_fill)
                    if mc_counts:
                        parts.append(mc_fn(mc_vals) if mc_vals is not None else mc_fill)
                    body.write("\t".join(parts) + "\n")
                    n_group += 1
                if n_group:
                    labels.append(label)
                    bounds.append(bounds[-1] + n_group)
        header["group_labels"] = labels
        header["group_boundaries"] = bounds
        if bounds[-1] == 0:
            raise SystemExit("matrix_join_all: no region has data in any column")
        with open(out_path, "wb") as out:
            out.write(gzip.compress(("@" + json.dumps(header, separators=(",", ":")) + "\n").encode(),
                                    compresslevel=COMPRESSLEVEL))
            with open(body_path, "rb") as body:
                while chunk := body.read(1 << 24):
                    out.write(chunk)
    finally:
        if os.path.exists(body_path):
            os.remove(body_path)
    return bounds[-1], dropped


def main(argv=None):
    ap = argparse.ArgumentParser(description="Join the most and mC matrices into the all matrix.")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--most", help="the most env's final matrix (non-mC)")
    ap.add_argument("--mc", help="the mC env's final matrix")
    args = ap.parse_args(argv)
    written, dropped = join_files(args.most, args.mc, args.output)
    print(f"matrix_join_all: {written} regions written, {dropped} with no data in any column dropped")


if __name__ == "__main__":
    main()
