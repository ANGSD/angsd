#!/bin/bash
# -doCounts test on four bams from smallBam, region 1:14000000-14002000.
# Each case writes the -dumpCounts files (and -doDepth/-doQsDist files), they are checked
# against the md5 sums in counts/counts.md5. The result must not depend on the number of threads.
# test5 (-dumpCounts 5 with more than one sample) used to crash, test6 checks that
# -setMinDepthInd also works with -minQfile (counts/minq.txt has a threshold per individual
# and base). -dumpCounts without -doCounts must give exit code 1.
if [[ ! $# -eq 2 ]]
then
	echo "Must supply an angsd binary and the folder with smallBam; will exit"
	exit 1
fi
ANGSD=${1}
BDIR=${2}

LOG=${0}.log
rm -f ${LOG}

MD5=md5sum
if [[ "$OSTYPE" == "darwin"* ]]; then
	MD5="$(cd "$(dirname "$0")" && pwd)/md5osx.sh"
fi

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD}" >> ${LOG}

TDIR=counts
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

for i in NA06985 NA06994 NA07000 NA07056;do
	ls ${BDIR}/small${i}.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam
done > ${RESDIR}/bam.filelist
BASE="-bam ${RESDIR}/bam.filelist -r 1:14000000-14002000 -doCounts 1"

RVAL=0
fail(){
	echo "--------------"
	echo "Problem with counts test: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

ARGS=("-dumpCounts 4"
      "-dumpCounts 2 -setMinDepthInd 2 -setMaxDepthInd 10"
      "-dumpCounts 3 -minInd 3 -setMinDepth 5 -setMaxDepth 40"
      "-dumpCounts 1 -doDepth 1 -maxDepth 20 -doQsDist 1 -minQ 0"
      "-dumpCounts 4 -setMaxDiffObs 1"
      "-dumpCounts 5"
      "-dumpCounts 4 -minQfile ${TDIR}/minq.txt -setMinDepthInd 2")

for n in 1 3; do
	D=${RESDIR}/threads${n}
	mkdir -p ${D}
	for t in ${!ARGS[@]}; do
		if ! ${ANGSD} ${BASE} ${ARGS[$t]} -nThreads ${n} -out ${D}/test${t} >> ${LOG} 2>&1;then
			fail "test${t} '${ARGS[$t]}' -nThreads ${n} exit code"
		fi
		for f in pos counts;do
			if [ -f ${D}/test${t}.${f}.gz ];then
				gunzip -f ${D}/test${t}.${f}.gz
			fi
		done
	done
	#counts.md5 has the md5 of every output file, test0.pos test0.counts ...
	if ! (cd ${D} && ${MD5} -c ../../counts.md5) >> ${LOG} 2>&1;then
		fail "output with -nThreads ${n} differs from counts.md5"
	fi
done

#options that need -doCounts
if ${ANGSD} -bam ${RESDIR}/bam.filelist -r 1:14000000-14002000 -dumpCounts 1 -out ${RESDIR}/bad >> ${LOG} 2>&1;then
	fail "-dumpCounts without -doCounts should give exit code 1"
fi

echo "==========$0===============" >> ${LOG}

exit $RVAL
