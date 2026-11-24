#!/bin/bash
# -beagle-GL reads a beagle file with genotype likelihoods (like -doGlf 2) as likelihoods,
# -beagle reads it as genotype probabilities. With the beagle file made from bamlikes/
# -doMaf 1, -doPost 1 and -doGeno 2 must give the same as with the bam files.
#
# ./testBeagleGL.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

RESDIR=beaglegl/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

ARGS=(-doMaf 1 -doPost 1 -doGeno 2)

${ANGSD} -bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doGlf 2 "${ARGS[@]}" -out ${RESDIR}/bam >> ${LOG} 2>&1 || { echo "angsd with bam input failed"; cat ${LOG}; exit 1; }
${ANGSD} -beagle-GL ${RESDIR}/bam.beagle.gz -fai thetacat/ref.fa.fai "${ARGS[@]}" -out ${RESDIR}/beagle >> ${LOG} 2>&1 || { echo "angsd with -beagle-GL input failed"; cat ${LOG}; exit 1; }

if [[ $(gzip -dc ${RESDIR}/bam.mafs.gz | wc -l) -lt 100 ]]; then
	echo "Problem with -beagle-GL: too few sites in ${RESDIR}/bam.mafs.gz"
	exit 1
fi
for f in mafs geno; do
	if ! diff <(gzip -dc ${RESDIR}/bam.${f}.gz) <(gzip -dc ${RESDIR}/beagle.${f}.gz) >> ${LOG}; then
		echo "Problem with -beagle-GL: ${f} differs from the bam input, see ${LOG}"
		exit 1
	fi
done

#the snp p-value must be a number (not nan) when a genotype likelihood in the file is 0
${ANGSD} -beagle-GL ${RESDIR}/bam.beagle.gz -fai thetacat/ref.fa.fai -doMaf 1 -SNP_pval 1e-2 -out ${RESDIR}/snp >> ${LOG} 2>&1 || { echo "angsd with -beagle-GL and -SNP_pval failed"; cat ${LOG}; exit 1; }
if [[ $(gzip -dc ${RESDIR}/snp.mafs.gz | tail -n +2 | wc -l) -eq 0 ]] || gzip -dc ${RESDIR}/snp.mafs.gz | cut -f6 | grep -q nan; then
	echo "Problem with -beagle-GL and -SNP_pval: no sites or a nan p-value"
	exit 1
fi
exit 0
