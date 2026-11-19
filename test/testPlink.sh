#!/bin/bash
# plink output test (-doPlink 2, tped/tfam)
# test1.beagle: 4 individuals, 7 sites on chromosomes 1 and 2, run with
# -postCutoff 0.5 -minInd 2. It has calls of all three genotypes, missing
# calls (written as 0 0), a call exactly at the cutoff (missing), a site
# with no calls (1_300) and a site with one call (1_400), both removed.
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

TDIR=plink
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

#the result must not depend on how many sites are read at a time
for n in 50 3; do
	${ANGSD} -beagle ${TDIR}/test1.beagle -fai hg19.fa.fai -doGeno -4 -doPlink 2 -postCutoff 0.5 -minInd 2 -nLines ${n} -out ${RESDIR}/test1_${n} >> ${LOG} 2>&1
	for ext in tped tfam; do
		if ! cmp ${RESDIR}/test1_${n}.${ext} ${TDIR}/test1.${ext} >> ${LOG} 2>&1;then
			echo "--------------"
			echo "Problem with plink test1 .${ext} -nLines ${n}"
			echo "--------------"
			cat ${LOG}
			echo "--------------"
			RVAL=1
		fi
	done
done

#invalid options must fail with a non-zero exit code
BADARGS=("-doGeno -4 -doPlink 1" "-doGeno -4 -doPlink 3" "-doPlink 2")
for args in "${BADARGS[@]}"; do
	${ANGSD} -beagle ${TDIR}/test1.beagle -fai hg19.fa.fai ${args} -out ${RESDIR}/bad >> ${LOG} 2>&1
	if [[ $? -eq 0 ]];then
		echo "--------------"
		echo "Problem with plink test: '${args}' did not fail"
		echo "--------------"
		cat ${LOG}
		echo "--------------"
		RVAL=1
	fi
done

echo "==========$0===============" >> ${LOG}

exit $RVAL
