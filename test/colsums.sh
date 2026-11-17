#!/bin/bash
#
# Print a short summary of a whitespace-separated table: the number of
# lines and fields, and the sum of every column from the second one on.
# Used to check large outputs (e.g. realSFS fst print) against a small
# reference with numcmp.sh, instead of storing the whole file.
#
#   ./colsums.sh FILE

if [ $# -ne 1 ] || [ ! -f "$1" ]; then
    echo "Usage: $0 FILE" >&2
    exit 2
fi

awk '
{
    if (NF > nf) nf = NF
    for (i = 2; i <= NF; i++) s[i] += $i
}
END {
    print "lines", NR
    print "fields", nf
    for (i = 2; i <= nf; i++) printf "sum%d %.6f\n", i, s[i]
}' "$1"
