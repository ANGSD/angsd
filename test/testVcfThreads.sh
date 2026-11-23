#!/bin/bash
# -doSaf on vcf input with threads (issue 598). The modules must not change
# chromosome before the threads are done with the chunks of the previous
# chromosome, otherwise the last sites of a chromosome are counted in the
# .saf.idx entry of the next. The saf with -nThreads 2 and 4 must be the same
# as with -nThreads 1, and every chromosome must have all its sites in order.
# vcfthreads/in.vcf.gz has 3 chromosomes with 500 sites each, -nLines 10
# gives many chunks per chromosome.
#
# ./testVcfThreads.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}
REALSFS=$(dirname ${ANGSD})/misc/realSFS

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD}" >> ${LOG}

if [[ ! -x ${REALSFS} ]]; then
	echo "Problem with -doSaf on vcf with threads: no ${REALSFS} to read the saf with"
	exit 1
fi

RESDIR=vcfthreads/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with -doSaf on vcf with threads: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

ARGS="-vcf-pl vcfthreads/in.vcf.gz -anc vcfthreads/ref.fa -doSaf 1 -doMajorMinor 1 -nLines 10"

for p in 1 2 4; do
	if ! ${ANGSD} ${ARGS} -nThreads ${p} -out ${RESDIR}/p${p} >> ${LOG} 2>&1; then
		problem "-nThreads ${p} exited with an error"
		continue
	fi
	${REALSFS} print ${RESDIR}/p${p}.saf.idx > ${RESDIR}/p${p}.txt 2>> ${LOG}
done

#the vcf has its 1500 sites sorted, so the saf positions must be the same in the same order
gunzip -c vcfthreads/in.vcf.gz | grep -v '^#' | cut -f1,2 > ${RESDIR}/vcfpos
if ! cut -f1,2 ${RESDIR}/p1.txt | cmp -s - ${RESDIR}/vcfpos; then
	problem "-nThreads 1: the saf sites are not the sites in the vcf"
fi
for p in 2 4; do
	if ! cmp -s ${RESDIR}/p1.txt ${RESDIR}/p${p}.txt; then
		problem "-nThreads ${p} gives a different saf than -nThreads 1"
		cut -f1 ${RESDIR}/p${p}.txt | uniq -c >> ${LOG}
	fi
done

echo "==========${0}===============" >> ${LOG}

exit $RVAL
