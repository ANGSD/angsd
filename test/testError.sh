#!/bin/bash
# -doError (#603).
# 1) The two samples in bamlikes/ with -minSites larger than the number of sites:
#    the sites left at the end were never used, so nothing was estimated. There
#    must now be one estimate.
# 2) With BAMDIR, 1:14000000-14030000 of 8 of the smallBam files in blocks of
#    500 sites: the threads added the sites to the same list at the same time,
#    which could crash. -nThreads 4 must give the same estimates as -nThreads 1,
#    and the sites after the last block must also be estimated.
# 3) With BAMDIR, speed and regression test: one estimate on the 10036 sites in
#    1:14000000-14070000 of all 33 smallBam files. The estimates must match
#    error/ref/speed.error.chunkunordered (made with this test), and the cpu time
#    is written to the log. It took about 32 s (Apple M-series), the test only
#    fails above MAXCPU seconds (default 300) so that slower machines pass.
#
# ./testError.sh angsd [BAMDIR]
# MAXCPU=120 ./testError.sh angsd BAMDIR     lower limit for the cpu time
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}
BAMDIR=${2}

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD} and BAMDIR: '${BAMDIR}'" >> ${LOG}

RESDIR=error/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with -doError: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

ARGS="-doMajorMinor 2 -doCounts 1 -doError 1"

#1)
if ! ${ANGSD} -bam bamlikes/bam.filelist ${ARGS} -minSites 100000 -out ${RESDIR}/few >> ${LOG} 2>&1; then
	problem "-minSites larger than the number of sites exited with an error"
elif [[ $(wc -l < ${RESDIR}/few.error.chunkunordered) -ne 1 ]]; then
	problem "-minSites larger than the number of sites didn't give one estimate"
fi

#wrong options must exit with 1, they used to exit with 0
for a in "-doMajorMinor 2 -doError 1" "-doMajorMinor 1 -doCounts 1 -doError 1"; do
	${ANGSD} -bam bamlikes/bam.filelist ${a} -out ${RESDIR}/bad >> ${LOG} 2>&1
	e=$?
	if [[ ${e} -ne 1 ]]; then
		problem "'${a}' exited with ${e}, not 1"
	fi
done

#2)
if [[ -n "${BAMDIR}" ]]; then
	ls ${BAMDIR}/*.bam | head -8 > ${RESDIR}/bam.filelist
	for p in 1 4; do
		${ANGSD} -bam ${RESDIR}/bam.filelist -r 1:14000000-14030000 ${ARGS} -minSites 500 -nThreads ${p} -out ${RESDIR}/t${p} > ${RESDIR}/t${p}.log 2>&1
		e=$?
		cat ${RESDIR}/t${p}.log >> ${LOG}
		if [[ ${e} -ne 0 ]]; then
			problem "-nThreads ${p} exited with ${e}"
		fi
	done
	n=$(grep -c "Triggering error estimation" ${RESDIR}/t1.log)
	if [[ ${n} -lt 2 ]]; then
		problem "less than 2 blocks of sites with -minSites 500 (${n})"
	fi
	if [[ $(wc -l < ${RESDIR}/t1.error.chunkunordered) -ne ${n} ]]; then
		problem "not one estimate for each block of sites"
	fi
	#all sites must be used, the ones after the last block of 500 included
	used=$(grep "Triggering error estimation" ${RESDIR}/t1.log | awk -F= '{s+=$2}END{print s}')
	kept=$(grep "Number of sites retained after filtering" ${RESDIR}/t1.log | awk '{print $NF}')
	if [[ "${used}" != "${kept}" ]]; then
		problem "${used} sites were used for the estimates, but ${kept} sites were kept"
	fi
	if ! cmp -s ${RESDIR}/t1.error.chunkunordered ${RESDIR}/t4.error.chunkunordered; then
		problem "-nThreads 4 doesn't give the same estimates as -nThreads 1"
	fi

	#3)
	ls ${BAMDIR}/*.bam > ${RESDIR}/all.filelist
	MAXCPU=${MAXCPU:-300}
	TIMEFORMAT=%U
	cpu=$( { time ${ANGSD} -bam ${RESDIR}/all.filelist -r 1:14000000-14070000 ${ARGS} -minSites 1000000 -out ${RESDIR}/speed > ${RESDIR}/speed.log 2>&1; } 2>&1 )
	e=$?
	cat ${RESDIR}/speed.log >> ${LOG}
	echo "-doError on 10036 sites x 33 individuals took ${cpu} s cpu (limit ${MAXCPU} s)" | tee -a ${LOG}
	if [[ ${e} -ne 0 ]]; then
		problem "the speed test exited with ${e}"
	elif ! grep -q "nsites=10036" ${RESDIR}/speed.log; then
		problem "the speed test didn't make one estimate on 10036 sites"
	elif ! ./numcmp.sh error/ref/speed.error.chunkunordered ${RESDIR}/speed.error.chunkunordered >> ${LOG} 2>&1; then
		problem "the estimates differ from error/ref/speed.error.chunkunordered"
	fi
	if ! awk -v c=${cpu} -v m=${MAXCPU} 'BEGIN{exit !(c<=m)}'; then
		problem "the speed test took ${cpu} s cpu, more than ${MAXCPU} s"
	fi
fi

exit ${RVAL}
