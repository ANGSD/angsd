#!/bin/bash
# -doGeno, -doPost and -doMaf stopped with exit code 0 when arguments were missing or wrong:
# -doGeno without -doPost, -geno_minDepth/-geno_maxDepth/-geno_minMM without -doCounts,
# -doPost 1 without -doMaf, -doPost without -GL, an unknown -doPost (after reading the data),
# -minMaf without -doMaf and a negative -SNP_pval. They must exit with 1.
# The same for -doMajorMinor: 1 without -GL, 2 without -doCounts, 3 without -sites, 4/7 without
# -ref and 5 without -anc, and an unknown -doMajorMinor and 7 without -ref crashed (segfault).
# An individual without reads has the same posterior for the 3 genotypes and must not be
# called, also with -postCutoff below 1/3 (it was called as major/major).
# -GL and -doGlf errors also exited with 0.
#
# ./testGenoArgs.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

RESDIR=genoargs/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

BAM="-bam bamlikes/bam.filelist"
IN="${BAM} -GL 1 -doMajorMinor 1"
RVAL=0
check1(){
	${ANGSD} "$@" -out ${RESDIR}/out >> ${LOG} 2>&1
	e=$?
	if [[ ${e} -ne 1 ]]; then
		echo "Problem with '$*': exit code ${e}, not 1"
		RVAL=1
	fi
}
check1 ${IN} -doMaf 1 -doGeno 2
check1 ${IN} -doMaf 1 -doPost 1 -doGeno 2 -geno_minDepth 2
check1 ${IN} -doMaf 1 -doPost 1 -doGeno 2 -geno_maxDepth 9
check1 ${IN} -doMaf 1 -doPost 1 -doGeno 2 -geno_minMM 0.9
check1 ${IN} -doPost 1 -doGeno 2
check1 ${BAM} -doMajorMinor 2 -doCounts 1 -doMaf 1 -doPost 1
check1 ${IN} -doMaf 1 -doPost 5
check1 ${IN} -doPost 2 -minMaf 0.05
check1 ${IN} -doMaf 1 -SNP_pval -0.5
check1 ${BAM} -doMajorMinor 1 -doMaf 1
check1 ${BAM} -GL 1 -doMajorMinor 2 -doMaf 1
for m in 3 4 5 7 9 -1; do
	check1 ${BAM} -GL 1 -doMajorMinor ${m} -doMaf 1
done
#-GL/-doGlf: unknown -GL, -doGlf without -GL, -GL 4/6 without -doCounts, -doGlf 2 without
#-doMajorMinor, a missing -errors file, the removed -doGlf 5 and an unknown -doGlf
check1 ${BAM} -GL 9
check1 ${BAM} -doGlf 1
check1 ${BAM} -GL 4
check1 ${BAM} -GL 6
check1 ${BAM} -GL 1 -doGlf 2
check1 ${BAM} -GL 4 -doCounts 1 -errors ${RESDIR}/notthere.errors
check1 ${BAM} -GL 1 -doGlf 5
check1 ${BAM} -GL 1 -doGlf 7

#and with everything needed it must work
${ANGSD} ${IN} -doMaf 1 -doPost 1 -doGeno 2 -doCounts 1 -geno_minDepth 2 -out ${RESDIR}/ok >> ${LOG} 2>&1
e=$?
if [[ ${e} -ne 0 ]]; then
	echo "Problem with -doGeno with all arguments: exit code ${e}, not 0"
	RVAL=1
fi
${ANGSD} ${BAM} -GL 1 -doMajorMinor 7 -ref thetacat/ref.fa -doMaf 1 -out ${RESDIR}/mm7 >> ${LOG} 2>&1
e=$?
if [[ ${e} -ne 0 ]]; then
	echo "Problem with -doMajorMinor 7 -ref: exit code ${e}, not 0"
	RVAL=1
fi

#individuals without reads must be -1 with -postCutoff 0.2, and the others must be called
${ANGSD} ${IN} -doMaf 1 -doPost 1 -doGeno 2 -postCutoff 0.2 -doCounts 1 -dumpCounts 2 -out ${RESDIR}/nodata >> ${LOG} 2>&1 || { echo "Problem with -postCutoff 0.2: exit code $?"; exit 1; }
gzip -dc ${RESDIR}/nodata.pos.gz | tail -n +2 | cut -f1,2 > ${RESDIR}/pos.txt
gzip -dc ${RESDIR}/nodata.counts.gz | tail -n +2 > ${RESDIR}/depth.txt
paste ${RESDIR}/pos.txt ${RESDIR}/depth.txt > ${RESDIR}/posdepth.txt
gzip -dc ${RESDIR}/nodata.geno.gz > ${RESDIR}/geno.txt
#chr pos depth1 depth2 ; chr pos geno1 geno2
res=$(awk 'NR==FNR{k=$1":"$2; d1[k]=$3; d2[k]=$4; next}
	{k=$1":"$2; if(!(k in d1)) next; n++
	 if(d1[k]==0) {z++; if($3!=-1) bad++} else if($3==-1) bad++
	 if(d2[k]==0) {z++; if($4!=-1) bad++} else if($4==-1) bad++}
	END{print n+0, z+0, bad+0}' ${RESDIR}/posdepth.txt ${RESDIR}/geno.txt)
read nsites nzero nbad <<< "${res}"
if [[ ${nsites} -lt 100 || ${nzero} -eq 0 ]]; then
	echo "Problem with the -postCutoff 0.2 test: ${nsites} sites, ${nzero} genotypes without reads"
	RVAL=1
elif [[ ${nbad} -ne 0 ]]; then
	echo "Problem with -postCutoff 0.2: ${nbad} genotypes called without reads or not called with reads"
	RVAL=1
fi
exit ${RVAL}
