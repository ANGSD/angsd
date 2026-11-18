#!/bin/bash
# bgen file reader test
# test1.bgen: 4 samples, 18 biallelic SNPs on chromosomes 1, 2 and 3,
# zlib compressed, 16 bits per probability, with sample identifiers.
if [[ ! $# -eq 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
else
	ANGSD=${1}
fi

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD}" >> ${LOG}

TDIR=bgen
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

#the result must not depend on how many sites are read at a time
for n in 50 3; do
	${ANGSD} -doMaf 4 -bgen ${TDIR}/test1.bgen -fai hg19.fa.fai -nLines ${n} -out ${RESDIR}/test1_${n} >> ${LOG} 2>&1
	if ! cmp <(gunzip -c ${RESDIR}/test1_${n}.mafs.gz) ${TDIR}/test1.mafs >> ${LOG} 2>&1;then
		echo "--------------"
		echo "Problem with bgen reader test1 -nLines ${n}"
		echo "--------------"
		cat ${LOG}
		echo "--------------"
		RVAL=1
	fi
done

echo "==========$0===============" >> ${LOG}

exit $RVAL
