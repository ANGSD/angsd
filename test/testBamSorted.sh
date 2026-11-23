#!/bin/bash
# Reads in files that are not sorted by position (samtools sort) were put on the
# wrong chromosome or dropped without a message, or angsd stopped with
# 'ReferenceID ... is not' or crashed, all with exit code 0 (#377).
# angsd must now stop with exit code 1 and say the file is not sorted.
# Makes two samples from bamlikes/s1.sam and s2.sam with the same reads on
# chr 1, 2 and 3 of thetacat/ref.fa, one sorted and one in different unsorted
# orders. Sorted files, also with unmapped reads at the end, must still work.
#
# ./testBamSorted.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD}" >> ${LOG}

RESDIR=bamsorted/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with unsorted bam files: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

hdr(){
	printf "@HD\tVN:1.0\tSO:coordinate\n"
	for c in 1 2 3; do printf "@SQ\tSN:${c}\tLN:150\n"; done
}
#the reads of sample $1 on chromosome $2
body(){
	grep -v '^@' bamlikes/$1.sam | awk -v c=$2 'BEGIN{OFS="\t"}{$3=c;print}'
}

{ hdr; body s2 1; body s2 2; body s2 3; } > ${RESDIR}/sorted.sam

#name and the reads of the second sample in that order
UNSORTED=("chr1_after_chr2"    "body s1 1 | head -300; body s1 2; body s1 1 | tail -300; body s1 3"
	  "chr1_after_chr3"    "body s1 1; body s1 2; body s1 3; body s1 1 | tail -5"
	  "chr2_inside_chr3"   "body s1 1; body s1 3 | head -300; body s1 2; body s1 3 | tail -300"
	  "within_chr1"        "body s1 1 | awk 'NR%2'; body s1 1 | awk 'NR%2==0'; body s1 2; body s1 3")

ARGS="-doCounts 1 -dumpCounts 2"
for ((i=0;i<${#UNSORTED[@]};i+=2)); do
	name=${UNSORTED[i]}
	{ hdr; eval "${UNSORTED[i+1]}"; } > ${RESDIR}/${name}.sam
	printf "${RESDIR}/sorted.sam\n${RESDIR}/${name}.sam\n" > ${RESDIR}/${name}.list
	${ANGSD} -b ${RESDIR}/${name}.list ${ARGS} -out ${RESDIR}/${name} > ${RESDIR}/${name}.log 2>&1
	e=$?
	cat ${RESDIR}/${name}.log >> ${LOG}
	if [[ ${e} -ne 1 ]]; then
		problem "unsorted file (${name}) exited with ${e}, not 1"
	fi
	if ! grep -q "is not sorted by position" ${RESDIR}/${name}.log; then
		problem "unsorted file (${name}) doesn't say the file is not sorted"
	fi

	#the same reads sorted must work
	{ hdr; grep -v '^@' ${RESDIR}/${name}.sam | sort -k3,3n -k4,4n; } > ${RESDIR}/${name}.sorted.sam
	printf "${RESDIR}/sorted.sam\n${RESDIR}/${name}.sorted.sam\n" > ${RESDIR}/${name}.sorted.list
	if ! ${ANGSD} -b ${RESDIR}/${name}.sorted.list ${ARGS} -out ${RESDIR}/${name}.sorted >> ${LOG} 2>&1; then
		problem "the sorted reads of ${name} exited with an error"
	fi
done

#unmapped reads without a position are at the end of a sorted file
{ cat ${RESDIR}/sorted.sam; grep -v '^@' bamlikes/s2.sam | head -3 | awk 'BEGIN{OFS="\t"}{$2=4;$3="*";$4=0;$6="*";print}'; } > ${RESDIR}/unmapped.sam
echo ${RESDIR}/unmapped.sam > ${RESDIR}/unmapped.list
if ! ${ANGSD} -b ${RESDIR}/unmapped.list ${ARGS} -out ${RESDIR}/unmapped >> ${LOG} 2>&1; then
	problem "a sorted file with unmapped reads at the end exited with an error"
fi

exit ${RVAL}
