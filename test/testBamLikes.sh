#!/bin/bash
# -GL 1 (bam_likes) test on two small samples in bamlikes/.
# s1.sam has sites with depth above 255, where angsd picks 255 bases at
# random with -seed, s2.sam is shallow. The result must not depend on
# the number of threads. Expected outputs are -doGlf 4 text.
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

TDIR=bamlikes
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

ARGS=("-seed 1" "-seed 1 -minQ 20 -minMapQ 20 -trim 3")
for t in 0 1; do
	for n in 1 3; do
		${ANGSD} -bam ${TDIR}/bam.filelist -GL 1 -doGlf 4 ${ARGS[$t]} -nThreads ${n} -out ${RESDIR}/test${t}_${n} >> ${LOG} 2>&1
		gunzip -c ${RESDIR}/test${t}_${n}.glf.gz > ${RESDIR}/test${t}_${n}.glf
		if ! ./numcmp.sh ${TDIR}/test${t}.glf ${RESDIR}/test${t}_${n}.glf >> ${LOG} 2>&1;then
			echo "--------------"
			echo "Problem with -GL 1 test${t} '${ARGS[$t]}' -nThreads ${n}"
			echo "--------------"
			cat ${LOG}
			echo "--------------"
			RVAL=1
		fi
	done
done

echo "==========$0===============" >> ${LOG}

exit $RVAL
