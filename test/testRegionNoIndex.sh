#!/bin/bash
# -r/-rf need an index. Without one (bamlikes/ has sam files) angsd printed an error but
# exited with 0 and left empty output files, so it must exit with 1.
#
# ./testRegionNoIndex.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

RESDIR=regionnoindex/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}
echo "1:60-140" > ${RESDIR}/regions.txt

RVAL=0
for a in "-r 1:60-140" "-rf ${RESDIR}/regions.txt" "-r 1:60-140 -nThreads 2"; do
	${ANGSD} -bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 ${a} -out ${RESDIR}/out >> ${LOG} 2>&1
	e=$?
	if [[ ${e} -ne 1 ]]; then
		echo "Problem with '${a}' on files without an index: exit code ${e}, not 1"
		RVAL=1
	fi
done
exit ${RVAL}
