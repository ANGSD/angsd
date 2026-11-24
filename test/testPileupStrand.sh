#!/bin/bash
# -pileup input must keep the strand of the reads with the reference base (#358).
# '.' is the reference base on the forward strand and ',' on the reverse strand, but
# both were read as the same base, so all reads with the reference base looked like
# forward strand reads and -sb_pval removed most sites.
# pileupstrand/bamlikes.pileup is the pileup of bamlikes/ made with
#   samtools mpileup -f thetacat/ref.fa -Q 0 -q 0 -B bamlikes/s1.sam bamlikes/s2.sam
# The strand counts and strand bias of -doSnpStat must be the same as with the bam files.
#
# ./testPileupStrand.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

RESDIR=pileupstrand/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

ARGS=(-doMajorMinor 2 -doCounts 1 -doSnpStat 1 -minQ 0)

${ANGSD} -bam bamlikes/bam.filelist -ref thetacat/ref.fa "${ARGS[@]}" -out ${RESDIR}/bam >> ${LOG} 2>&1 || { echo "angsd with bam input failed"; cat ${LOG}; exit 1; }
${ANGSD} -pileup pileupstrand/bamlikes.pileup -fai thetacat/ref.fa.fai -nInd 2 "${ARGS[@]}" -out ${RESDIR}/pileup >> ${LOG} 2>&1 || { echo "angsd with pileup input failed"; cat ${LOG}; exit 1; }

#chr pos, the counts +Major +Minor -Major -Minor and the strand bias. The other columns
#use the mapping quality and the position in the read, which a pileup doesn't have
gzip -dc ${RESDIR}/bam.snpStat.gz | cut -f1-4 > ${RESDIR}/bam.strand
gzip -dc ${RESDIR}/pileup.snpStat.gz | cut -f1-4 > ${RESDIR}/pileup.strand
if [[ $(wc -l < ${RESDIR}/bam.strand) -lt 100 ]]; then
	echo "Problem with -pileup strand: too few sites in ${RESDIR}/bam.strand"
	exit 1
fi
if ! diff ${RESDIR}/bam.strand ${RESDIR}/pileup.strand >> ${LOG}; then
	echo "Problem with -pileup strand: the strand counts differ from the bam input, see ${LOG}"
	exit 1
fi
exit 0
