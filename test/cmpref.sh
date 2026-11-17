#!/bin/bash
#
# Compare test output against the reference files in a directory, with
# the tolerance of numcmp.sh. What is compared depends on the name of
# each reference file:
#
#   NAME.sums   the column sums of OUTDIR/NAME (see colsums.sh)
#   NAME        the whole file OUTDIR/NAME
#
#   ./cmpref.sh REFDIR OUTDIR [REL_TOL]
#
# REL_TOL is the relative tolerance for whole files (default: that of
# numcmp.sh). The .sums comparisons always use their own tolerance.
# Exits 1 if any comparison fails, after checking all of them.

if [ $# -lt 2 ] || [ $# -gt 3 ] || [ ! -d "$1" ]; then
    echo "Usage: $0 REFDIR OUTDIR [REL_TOL]" >&2
    exit 2
fi

here=$(cd "$(dirname "$0")" && pwd)
refdir=$1
outdir=$2
reltol=$3
tmp=$(mktemp -d "${TMPDIR:-/tmp}/cmpref.XXXXXX") || exit 2
trap 'rm -rf "$tmp"' EXIT

status=0
for ref in "$refdir"/*; do
    [ -f "$ref" ] || continue
    name=$(basename "$ref")
    ok=1
    case "$name" in
	*.sums)
	    out="$outdir/${name%.sums}"
	    [ -f "$out" ] || { echo "cmpref: missing file '$out'"; status=1; continue; }
	    "$here/colsums.sh" "$out" > "$tmp/$name"
	    #sums add up the last-digit noise of many values, so allow more
	    "$here/numcmp.sh" "$ref" "$tmp/$name" 1e-4 || ok=0
	    ;;
	*)
	    if [ -n "$reltol" ]; then
		"$here/numcmp.sh" "$ref" "$outdir/$name" 1e-5 "$reltol" || ok=0
	    else
		"$here/numcmp.sh" "$ref" "$outdir/$name" || ok=0
	    fi
	    ;;
    esac
    if [ $ok -eq 1 ]; then echo "cmpref: $name OK"; else status=1; fi
done
exit $status
