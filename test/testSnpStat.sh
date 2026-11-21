#!/bin/bash
# -doSnpStat test on 12 of the smallBam files (needs BAMDIR). The SNPs in
# 1:14000000-14030000 are run with each filter, and the .snpStat and .mafs
# output must match snpstat/expected/, which was made with
# angsd 0.943-10-g1c624c2.
#
# ./testSnpStat.sh angsd BAMDIR             runs the test
# ./testSnpStat.sh angsd BAMDIR expected    remakes snpstat/expected/ with this angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}
BAMDIR=${2}
MAKEEXP=0
if [[ $# -eq 3 && "$3" == "expected" ]]; then
	MAKEEXP=1
fi

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD} and BAMDIR: '${BAMDIR}'" >> ${LOG}

if [[ -z "${BAMDIR}" ]]; then
	echo "No BAMDIR, -doSnpStat is not tested" >> ${LOG}
	exit 0
fi

TDIR=snpstat
EXPDIR=${TDIR}/expected
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with -doSnpStat: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

ls ${BAMDIR}/*.bam | sort | head -12 > ${RESDIR}/bam.filelist
BASE="-bam ${RESDIR}/bam.filelist -r 1:14000000-14030000 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doCounts 1 -doSnpStat 1"
NAMES=(hwe sb edge hwefilter qscore mapq hetbias)
ARGS=("-doHWE 1"
      "-doHWE 1 -sb_pval 0.05"
      "-doHWE 1 -edge_pval 0.5"
      "-doHWE 1 -hwe_pval 0.05"
      "-doHWE 1 -qscore_pval 0.05"
      "-doHWE 1 -mapQ_pval 0.05"
      "-doHWE 1 -doPost 1 -doGeno 1 -hetbias_pval 0.05")
for i in ${!NAMES[@]}; do
	n=${NAMES[$i]}
	if ! ${ANGSD} ${BASE} ${ARGS[$i]} -out ${RESDIR}/${n} >> ${LOG} 2>&1; then
		problem "'${ARGS[$i]}' exited with an error"
		continue
	fi
	for t in snpStat mafs; do
		gunzip -c ${RESDIR}/${n}.${t}.gz > ${RESDIR}/${n}.${t}
		if [[ ${MAKEEXP} -eq 1 ]]; then
			cp ${RESDIR}/${n}.${t} ${EXPDIR}/${n}.${t}
		elif ! ./numcmp.sh ${EXPDIR}/${n}.${t} ${RESDIR}/${n}.${t} >> ${LOG} 2>&1; then
			problem "the .${t} for '${ARGS[$i]}' differs from ${EXPDIR}/${n}.${t}"
		fi
	done
done
if [[ ${MAKEEXP} -eq 1 ]]; then
	echo "Made the files in ${EXPDIR}"
fi

echo "==========$0===============" >> ${LOG}

exit $RVAL
