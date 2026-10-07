#!/usr/bin/env bash
# merge_cx_reports.sh — sum replicate CX reports into one.
#
# Usage: merge_cx_reports.sh <threads> <out.CX_report.txt.gz> <in.CX_report.txt.gz>...
#
# Output has the CX report's seven columns, with the methylated and
# unmethylated counts summed over the replicates at each cytosine.
#
# Bismark's genome-wide cytosine report lists every cytosine in the genome,
# covered or not, so replicate reports hold the same rows for each
# chromosome. (Not in the same chromosome order: Bismark 3.x writes
# chromosomes in a different order on every run.) Each report is split into
# per-chromosome blocks, and the blocks are summed side by side -- linear
# throughout. Reports whose rows don't line up (dmC reports carry only
# covered positions) are merged by sorting them together and summing with
# bedtools merge, which on a 1 Gb genome runs for many hours.

set -euo pipefail

if [[ $# -lt 3 ]]; then
    printf "Usage: %s <threads> <out.CX_report.txt.gz> <in.CX_report.txt.gz>...\n" "$0" >&2
    exit 1
fi

threads="$1"
out="$2"
shift 2
n=$#
split_awk="$(dirname "$0")/split_by_chrom.awk"

work=$(mktemp -d "${TMPDIR:-/tmp}/mergecx.XXXXXX")
trap 'rm -rf "$work"' EXIT

merge_by_sort() {
    printf "merge_cx_reports: replicate reports do not line up row for row; merging by sort\n"
    rm -rf "$work"/rep*
    zcat "$@" \
        | sort -k1,1 -k2,2n \
        | awk -v OFS='\t' '{print $1, $2-1, $2, $3, $4, $5, $6, $7}' > "$work/merged.bed"
    bedtools merge -d -1 -o distinct,sum,sum,distinct,distinct -c 4,5,6,7,8 -i "$work/merged.bed" \
        | awk -v OFS='\t' '{print $1, $3, $4, $5, $6, $7, $8}' \
        | pigz -p "$threads" > "$out"
}

# Split every report at once. A report that is not grouped by chromosome
# with ascending positions exits 3.
pids=()
for i in $(seq 1 "$n"); do
    mkdir "$work/rep$i"
    zcat "${!i}" | LC_ALL=C awk -v dir="$work/rep$i" -f "$split_awk" &
    pids+=($!)
done
grouped=0
for pid in "${pids[@]}"; do
    rc=0
    wait "$pid" || rc=$?
    if [[ "$rc" -eq 3 ]]; then
        grouped=3
    elif [[ "$rc" -ne 0 ]]; then
        exit "$rc"
    fi
done
if [[ "$grouped" -ne 0 ]]; then
    merge_by_sort "$@"
    exit 0
fi

# Each block file of replicate 1 alongside the same chromosome's block from
# every other replicate. A chromosome missing from any replicate means the
# reports don't hold the same rows.
blocks="$work/blocks"
if ! awk -F'\t' -v n="$n" '
        { path[FILENAME, $1] = $2 }
        FILENAME == ARGV[1] { order[++m] = $1 }
        END {
            for (r = 2; r <= n; r++) if (count(ARGV[r]) != m) exit 3
            for (i = 1; i <= m; i++) {
                line = ""
                for (r = 1; r <= n; r++) {
                    if (!((ARGV[r], order[i]) in path)) exit 3
                    line = line (r > 1 ? "\t" : "") path[ARGV[r], order[i]]
                }
                print line
            }
        }
        function count(f,   k, c, parts) {
            c = 0
            for (k in path) { split(k, parts, SUBSEP); if (parts[1] == f) c++ }
            return c
        }' "$work"/rep*/index > "$blocks"; then
    merge_by_sort "$@"
    exit 0
fi

# A row whose chromosome, position or strand differs between replicates
# exits 3; so does a block running out early, since paste then fills its
# columns with empty fields.
summed=0
while IFS=$'\t' read -r -a files; do
    paste "${files[@]}"
done < "$blocks" \
    | awk -F'\t' -v OFS='\t' -v n="$n" '
        {
            m = $4; u = $5
            for (k = 1; k < n; k++) {
                o = 7 * k
                if ($(o+1) != $1 || $(o+2) != $2 || $(o+3) != $3) exit 3
                m += $(o+4); u += $(o+5)
            }
            print $1, $2, $3, m, u, $6, $7
        }' \
    | pigz -p "$threads" > "$out" || summed=$?

if [[ "$summed" -ne 0 ]]; then
    merge_by_sort "$@"
fi
