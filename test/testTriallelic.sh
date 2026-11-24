#!/bin/bash
# -rmTriallelic must keep biallelic sites and remove triallelic sites.
# addProtectN left out the first genotype (AA) when summing the likelihoods, so biallelic
# sites with the allele A could be removed.
# triallelic/sites.glf is made by triallelic/mk.py: site 1 is biallelic A/C, site 2 biallelic
# G/T and site 3 triallelic C/G/T.
#
# ./testTriallelic.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

RESDIR=triallelic/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

${ANGSD} -glf10_text triallelic/sites.glf -nInd 10 -fai thetacat/ref.fa.fai -doMajorMinor 1 -doMaf 1 -rmTriallelic 0.05 -out ${RESDIR}/tri >> ${LOG} 2>&1 || { echo "angsd -rmTriallelic failed"; cat ${LOG}; exit 1; }
kept=$(gzip -dc ${RESDIR}/tri.mafs.gz | tail -n +2 | cut -f2 | tr '\n' ' ')
if [[ "${kept}" != "1 2 " ]]; then
	echo "Problem with -rmTriallelic: kept the sites '${kept}', expected '1 2 '"
	exit 1
fi
exit 0
