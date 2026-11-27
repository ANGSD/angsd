#!/bin/bash
# -show 1 with an insertion or deletion before the first aligned base of a read (5S2D10M,
# 5H2D10M, 2I2D10M). The indel was bound to the position before the read: it wrote a
# '-2NN' to a read that doesn't belong there, or outside the nodes when no read covered
# that position (segfault for first.sam, where such a read is the first one).
# showindel/*.samtools are made with
#   samtools mpileup -x -B -Q 0 -q 0 -f showindel/ref.fa showindel/<name>.sam
# and -show 1 -ref must give the same.
#
# ./testShowIndel.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

TDIR=showindel
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0
for f in first indels; do
	echo ${TDIR}/${f}.sam > ${RESDIR}/${f}.filelist
	${ANGSD} -b ${RESDIR}/${f}.filelist -show 1 -ref ${TDIR}/ref.fa -out ${RESDIR}/${f} > ${RESDIR}/${f}.show 2>> ${LOG}
	e=$?
	if [[ ${e} -ne 0 ]]; then
		echo "Problem with -show 1 on ${f}.sam: exit code ${e}"
		RVAL=1
	elif ! diff ${RESDIR}/${f}.show ${TDIR}/${f}.samtools >> ${LOG}; then
		echo "Problem with -show 1 on ${f}.sam: differs from samtools mpileup, see ${LOG}"
		RVAL=1
	fi
done
exit ${RVAL}
