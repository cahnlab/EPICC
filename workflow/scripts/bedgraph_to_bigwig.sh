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
# The input does not need to be sorted. Input that is already grouped by
# chromosome with ascending starts -- every bedGraph the pipeline builds, since
# they come from CX reports and coordinate-sorted BAMs -- is put in order by
# concatenating per-chromosome blocks, which is linear. Anything else goes
# through sort, which on a 10 GB CHH track from a 1 Gb genome took hours.

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

mkdir "$work/chr"

# Split the data into one file per chromosome; exit 3 means the input is not
# grouped by chromosome with ascending starts, and is sorted instead.
grouped=0
LC_ALL=C awk -v dir="$work/chr" -f "$(dirname "$0")/split_by_chrom.awk" "$in_bg" \
    || grouped=$?
if [[ "$grouped" -ne 0 && "$grouped" -ne 3 ]]; then
    exit "$grouped"
fi

# Chromosomes in chrom.sizes with no interval in the data, padded with a
# single zero-value base. Keyed on FILENAME rather than NR==FNR so an empty
# bedGraph (empty index) still pads every chromosome.
if [[ "$grouped" -eq 0 ]]; then
    seen_from="$work/chr/index"
else
    seen_from="$in_bg"
fi
awk -v first="$seen_from" -v OFS='\t' '
    FILENAME == first { seen[$1] = 1; next }
    !($1 in seen) { print $1, 0, 1, 0 }
' "$seen_from" "$chrom_sizes" > "$work/pad.bedGraph"

n_pad=$(wc -l < "$work/pad.bedGraph")
if [[ "$n_pad" -gt 0 ]]; then
    printf "bedgraph_to_bigwig: %s has no coverage on %d chromosome(s); padding with zeros\n" \
           "$(basename "$in_bg")" "$n_pad"
fi

if [[ "$grouped" -eq 0 ]]; then
    # Each pad line is a block of its own, so it sorts in with the rest.
    n=$(wc -l < "$work/chr/index")
    while IFS=$'\t' read -r chrom _ _ _; do
        n=$((n + 1))
        printf '%s\t0\t1\t0\n' "$chrom" > "$work/chr/$n"
        printf '%s\t%s\n' "$chrom" "$work/chr/$n" >> "$work/chr/index"
    done < "$work/pad.bedGraph"
    LC_ALL=C sort -t $'\t' -k1,1 "$work/chr/index" | cut -f2 | xargs -r -d '\n' cat \
        > "$work/sorted.bedGraph"
else
    cat "$in_bg" "$work/pad.bedGraph" \
        | LC_COLLATE=C sort -k1,1 -k2,2n > "$work/sorted.bedGraph"
fi
rm -rf "$work/chr"

bedGraphToBigWig "$work/sorted.bedGraph" "$chrom_sizes" "$out_bw"
