#!/bin/bash
# -doSaf with filters that run after abcSaf (here -doSnpStat with -hwe_pval,
# -sb_pval, -edge_pval, -qscore_pval and -mapQ_pval). The saf of a site does
# not depend on the other sites, so every site in the filtered .saf must have
# the same values as in a run without the filter, and the sites must be the
# ones in the .mafs of the filtered run.
# Runs on the two samples in bamlikes/ with -doSaf 1 (thetacat/ref.fa as -anc)
# and -doSaf 5.
#
# ./testSafFilter.sh angsd
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
	echo "Problem with -doSaf and filters: no ${REALSFS} to read the saf with"
	exit 1
fi

RESDIR=saffilter/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with -doSaf and filters: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

FILTERS=("-doHWE 1 -doSnpStat 1 -hwe_pval 0.5"
	 "-doSnpStat 1 -sb_pval 0.5"
	 "-doSnpStat 1 -edge_pval 0.5"
	 "-doSnpStat 1 -qscore_pval 0.5"
	 "-doSnpStat 1 -mapQ_pval 0.5")

for dosaf in 1 5; do
	BASE="-bam bamlikes/bam.filelist -GL 1 -doSaf ${dosaf} -doMajorMinor 1 -doMaf 1 -doCounts 1"
	if [[ ${dosaf} -eq 1 ]]; then
		BASE="${BASE} -anc thetacat/ref.fa"
	fi
	all=${RESDIR}/saf${dosaf}_all
	if ! ${ANGSD} ${BASE} -out ${all} >> ${LOG} 2>&1; then
		problem "-doSaf ${dosaf} without filters exited with an error"
		continue
	fi
	${REALSFS} print ${all}.saf.idx > ${all}.txt 2>> ${LOG}
	nall=$(wc -l < ${all}.txt)
	for i in ${!FILTERS[@]}; do
		f=${RESDIR}/saf${dosaf}_f${i}
		if ! ${ANGSD} ${BASE} ${FILTERS[$i]} -out ${f} >> ${LOG} 2>&1; then
			problem "-doSaf ${dosaf} '${FILTERS[$i]}' exited with an error"
			continue
		fi
		${REALSFS} print ${f}.saf.idx > ${f}.txt 2>> ${LOG}
		gunzip -c ${f}.mafs.gz | tail -n +2 | cut -f1,2 > ${f}.mafspos
		n=$(wc -l < ${f}.txt)
		echo "-doSaf ${dosaf} '${FILTERS[$i]}': ${n} of ${nall} sites" >> ${LOG}
		if [[ ${n} -eq 0 || ${n} -ge ${nall} ]]; then
			problem "-doSaf ${dosaf} '${FILTERS[$i]}' kept ${n} of ${nall} sites, the filter is not tested"
		fi
		if ! cut -f1,2 ${f}.txt | cmp -s - ${f}.mafspos; then
			problem "-doSaf ${dosaf} '${FILTERS[$i]}': the saf sites are not the sites in the .mafs"
		fi
		if ! awk -F'\t' 'NR==FNR{a[$1":"$2]=$0;next}{if(a[$1":"$2]!=$0){print "differs: "$1":"$2;bad++}}END{exit bad>0}' ${all}.txt ${f}.txt >> ${LOG}; then
			problem "-doSaf ${dosaf} '${FILTERS[$i]}': saf values differ from the run without the filter"
		fi
	done
done

echo "==========$0===============" >> ${LOG}

exit $RVAL
