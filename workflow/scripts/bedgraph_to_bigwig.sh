#!/usr/bin/env bash
# bedgraph_to_bigwig.sh — bedGraph → bigWig with a complete chromosome header.
#
# Usage: bedgraph_to_bigwig.sh <in.bedGraph> <chrom.sizes> <out.bw>
#
# bedGraphToBigWig only writes the chromosomes that appear in its input, so a
# track with no coverage on a chromosome produces a bigWig whose header omits
# it. deeptools 4 treats a query against a missing chromosome as a hard error
# ("The passed chromosome (Chr2) was incorrect") and aborts computeMatrix, so
# every chromosome absent from the data is padded here with a single
# zero-value base. An empty input yields an all-zero bigWig over the whole
# genome rather than a one-chromosome stub.
#
# The input does not need to be sorted; this sorts before conversion.

set -euo pipefail

if [[ $# -ne 3 ]]; then
    printf "Usage: %s <in.bedGraph> <chrom.sizes> <out.bw>\n" "$0" >&2
    exit 1
fi

in_bg="$1"
chrom_sizes="$2"
out_bw="$3"

if [[ ! -f "$in_bg" ]]; then
    printf "ERROR: bedGraph not found: %s\n" "$in_bg" >&2
    exit 1
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/bg2bw.XXXXXX")
trap 'rm -rf "$work"' EXIT

# Chromosomes in chrom.sizes with no interval in the data. Keyed on FILENAME
# rather than NR==FNR so an empty bedGraph still pads every chromosome.
awk -v first="$in_bg" -v OFS='\t' '
    FILENAME == first { seen[$1] = 1; next }
    !($1 in seen) { print $1, 0, 1, 0 }
' "$in_bg" "$chrom_sizes" > "$work/pad.bedGraph"

n_pad=$(wc -l < "$work/pad.bedGraph")
if [[ "$n_pad" -gt 0 ]]; then
    printf "bedgraph_to_bigwig: %s has no coverage on %d chromosome(s); padding with zeros\n" \
           "$(basename "$in_bg")" "$n_pad"
fi

cat "$in_bg" "$work/pad.bedGraph" \
    | LC_COLLATE=C sort -k1,1 -k2,2n > "$work/sorted.bedGraph"

bedGraphToBigWig "$work/sorted.bedGraph" "$chrom_sizes" "$out_bw"
