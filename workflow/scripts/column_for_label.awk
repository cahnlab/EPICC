# column_for_label.awk — print the 1-based index of the header column named
# <label> in a multiBigwigSummary --outRawCounts table.
#
# Usage: awk -v label=<label> -f column_for_label.awk <values.tab>
#
# deeptools quotes every header name ('WT_H3K9me2') and prefixes the first
# with '#', so names are unquoted before an exact comparison. The labels and
# the header both come from the same -l list, so a label matches exactly one
# column; any other count exits 1 rather than picking a column.

NR == 1 {
    for (i = 1; i <= NF; i++) {
        name = $i
        gsub(/^#?'|'$/, "", name)
        if (name == label) { col = i; n++ }
    }
    if (n != 1) {
        printf "column_for_label: '%s' matches %d columns in %s\n", label, n, FILENAME > "/dev/stderr"
        exit 1
    }
    print col
    exit
}
