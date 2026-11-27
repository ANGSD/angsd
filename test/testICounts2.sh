#!/bin/bash
# -iCounts 2 skips sites with a deletion or an insertion. The nodes are reused from a pool,
# and a reused node kept the deletions of its last site, so -iCounts 2 also skipped sites
# without indels (320 of 92962 sites on chr 1 of the file below).
# The sites of -iCounts 2 must be the sites of -iCounts 1 without the sites with a '*'
# (deletion) or a '+' (insertion after the base) in -show 1, which is the same as samtools.
# -iCounts other than 0, 1 and 2 must exit with 1.
#
# ./testICounts2.sh angsd BAMDIR
if [[ ! $# -eq 2 ]]
then
	echo "Must supply an angsd binary and the folder with smallBam; will exit"
	exit 1
fi
ANGSD=${1}
BDIR=${2}

LOG=${0}.log
rm -f ${LOG}

RESDIR=icounts2/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

ls ${BDIR}/smallNA06985.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam > ${RESDIR}/bam.filelist
ARGS="-b ${RESDIR}/bam.filelist -r 1 -minQ 0"

#positions of an .icnts.gz: each site is 5 ints, the position and the 4 counts
icntsPos(){
	gzip -dc $1 | od -An -v -t d4 | tr -s ' \t' '\n\n' | grep -v '^$' | awk 'NR%5==1' | sort
}

RVAL=0
${ANGSD} ${ARGS} -doCounts 1 -iCounts 1 -out ${RESDIR}/i1 >> ${LOG} 2>&1 || { echo "Problem with -iCounts 1: exit code $?"; exit 1; }
${ANGSD} ${ARGS} -doCounts 1 -iCounts 2 -out ${RESDIR}/i2 >> ${LOG} 2>&1 || { echo "Problem with -iCounts 2: exit code $?"; exit 1; }
${ANGSD} ${ARGS} -show 1 -out ${RESDIR}/show > ${RESDIR}/show.txt 2>> ${LOG} || { echo "Problem with -show 1: exit code $?"; exit 1; }

icntsPos ${RESDIR}/i1.icnts.gz > ${RESDIR}/all.pos
icntsPos ${RESDIR}/i2.icnts.gz > ${RESDIR}/i2.pos
awk '$4 ~ /\*/ || $4 ~ /\+[0-9]/ {print $2}' ${RESDIR}/show.txt | sort > ${RESDIR}/indel.pos
comm -23 ${RESDIR}/all.pos ${RESDIR}/indel.pos > ${RESDIR}/expected.pos

if [[ $(wc -l < ${RESDIR}/all.pos) -lt 1000 || $(wc -l < ${RESDIR}/indel.pos) -eq 0 ]]; then
	echo "Problem with -iCounts 2 test: too few sites or no indels"
	RVAL=1
elif ! cmp ${RESDIR}/expected.pos ${RESDIR}/i2.pos >> ${LOG} 2>&1; then
	echo "Problem with -iCounts 2: the sites are not the sites of -iCounts 1 without indels ($(comm -23 ${RESDIR}/expected.pos ${RESDIR}/i2.pos | wc -l) sites missing)"
	RVAL=1
fi

#other values than 0, 1 and 2 made an empty .icnts.gz and exit code 0, they must exit with 1
for v in 3 -1; do
	${ANGSD} ${ARGS} -doCounts 1 -iCounts ${v} -out ${RESDIR}/bad >> ${LOG} 2>&1
	e=$?
	if [[ ${e} -ne 1 ]]; then
		echo "Problem with -iCounts ${v}: exit code ${e}, not 1"
		RVAL=1
	fi
done
exit ${RVAL}
