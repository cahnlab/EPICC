# split_by_chrom.awk — split a chromosome/position-keyed file into one file
# per chromosome.
#
# Usage: LC_ALL=C awk -v dir=<dir> -f split_by_chrom.awk <file>
#
# Writes each chromosome's lines to <dir>/<n> and lists "<chrom>\t<dir>/<n>"
# in <dir>/index, in order of appearance. The input must already be grouped
# by chromosome with ascending positions (column 2) -- true of CX reports and
# of anything built from them -- and then ordering it is just concatenating
# the blocks, which is linear where sort is n log n. Input that is not
# grouped that way exits 3, for the caller to fall back to sort.

BEGIN { OFS = "\t"; printf "" > (dir "/index") }

$1 != chrom {
    if ($1 in seen) exit 3
    close(out)
    chrom = $1; seen[chrom] = 1; out = dir "/" ++n; last = -1
    print chrom, out > (dir "/index")
}

$2 + 0 < last { exit 3 }

{ last = $2 + 0; print > out }
