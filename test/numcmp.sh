#!/bin/bash
#
# Compare two whitespace-separated text files, allowing small numeric
# differences. Used instead of md5 for outputs whose last printed digit
# can differ between platforms (e.g. Linux and macOS math libraries).
#
#   ./numcmp.sh REFERENCE NEW [ABS_TOL] [REL_TOL]
#
# The files must have the same number of lines and fields per line, and
# every non-numeric field (names, headers, inf, nan) must match exactly.
# Two numbers match if they differ by at most ABS_TOL (default 1e-5) or
# by at most REL_TOL (default 1e-9) relative to the larger of the two.
# Prints up to 10 mismatches and exits 1 if there are any.

if [ $# -lt 2 ]; then
    echo "Usage: $0 REFERENCE NEW [ABS_TOL] [REL_TOL]" >&2
    exit 2
fi

ref=$1
new=$2
abstol=${3:-1e-5}
reltol=${4:-1e-9}

for f in "$ref" "$new"; do
    if [ ! -f "$f" ]; then
	echo "numcmp: missing file '$f'" >&2
	exit 1
    fi
done
if [ ! -s "$ref" ]; then
    echo "numcmp: reference '$ref' is empty" >&2
    exit 1
fi

awk -v abstol="$abstol" -v reltol="$reltol" -v ref="$ref" -v new="$new" '
function isnum(s) {
    return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/
}
function abs(x) { return x < 0 ? -x : x }
function report(msg) {
    bad++
    if (bad <= 10) print "numcmp: " new " line " FNR ": " msg
}
NR == FNR { line[FNR] = $0; nref = FNR; next }
{
    if (FNR > nref) { report("extra line"); next }
    n = split(line[FNR], r)
    if (n != NF) { report("has " NF " fields, reference has " n); next }
    for (i = 1; i <= NF; i++) {
	a = r[i]; b = $i
	if (a == b) continue
	if (isnum(a) && isnum(b)) {
	    d = abs(a - b); m = abs(a) > abs(b) ? abs(a) : abs(b)
	    if (d <= abstol + 0 || d <= (reltol + 0) * m) continue
	}
	report("field " i " is " b ", reference has " a)
    }
}
END {
    if (FNR < nref && NR > nref) { bad++; print "numcmp: " new " has " FNR " lines, reference has " nref }
    if (NR == nref && nref > 0) { bad++; print "numcmp: " new " is empty" }
    if (bad > 10) print "numcmp: " bad " mismatches in total"
    exit (bad > 0)
}' "$ref" "$new"
