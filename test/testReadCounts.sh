#!/bin/bash
# The read line at the end of a run with BAM/CRAM input (issue 474). A paired
# read without the proper pair flag is removed by -only_proper_pairs 1 (the
# default), and angsd must say so instead of only reporting 0 sites. The used
# and removed reads must add up to the reads read.
#
# ./testReadCounts.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD}" >> ${LOG}

RESDIR=readcounts/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with the read counts: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

#one read with flag 161: paired, mate reverse, second in pair, not properly paired
SAM=${RESDIR}/notproper.sam
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n' > ${SAM}
printf 'r1\t161\tchr1\t100\t60\t30M\t=\t300\t230\tGGCATCGACATGTGCAGCGACCCTTGCGAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' >> ${SAM}
echo ${SAM} > ${RESDIR}/notproper.list

${ANGSD} -b ${RESDIR}/notproper.list -GL 1 -doMajorMinor 1 -doMaf 1 -out ${RESDIR}/proper1 > ${RESDIR}/proper1.log 2>&1
cat ${RESDIR}/proper1.log >> ${LOG}
if ! grep -q 'Total number of reads read: 1$' ${RESDIR}/proper1.log || ! grep -q 'Number of reads used after the read filters: 0$' ${RESDIR}/proper1.log || ! grep -q '1 removed by -only_proper_pairs$' ${RESDIR}/proper1.log; then
	problem "-only_proper_pairs 1: the read counts do not say the read was removed"
fi
if ! grep -q 'use -only_proper_pairs 0 to keep them' ${RESDIR}/proper1.log; then
	problem "-only_proper_pairs 1: no hint to use -only_proper_pairs 0"
fi
if ! grep -q 'Number of reads used after the read filters: 0$' ${RESDIR}/proper1.arg; then
	problem "-only_proper_pairs 1: no read counts in the .arg file"
fi

${ANGSD} -b ${RESDIR}/notproper.list -GL 1 -doMajorMinor 1 -doMaf 1 -only_proper_pairs 0 -out ${RESDIR}/proper0 > ${RESDIR}/proper0.log 2>&1
cat ${RESDIR}/proper0.log >> ${LOG}
if ! grep -q 'Number of reads used after the read filters: 1$' ${RESDIR}/proper0.log; then
	problem "-only_proper_pairs 0: the read is not used"
fi
if grep -q 'sites analyzed: 0$' ${RESDIR}/proper0.log; then
	problem "-only_proper_pairs 0: no sites"
fi

#the used and removed reads add up to the reads read
${ANGSD} -b bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -minMapQ 20 -out ${RESDIR}/bamlikes > ${RESDIR}/bamlikes.log 2>&1
cat ${RESDIR}/bamlikes.log >> ${LOG}
if ! awk '/Total number of reads read:/{nread=$NF}
	/Number of reads used after the read filters:/{used=$NF}
	/^\t\t *[0-9]+ /{removed+=$1; nfilt++}
	END{print "read "nread" used "used" removed "removed; exit !(nread>0 && nread==used+removed && nfilt>0)}' ${RESDIR}/bamlikes.log >> ${LOG}; then
	problem "bamlikes: the used and removed reads do not add up to the reads read, or no reads were removed by -minMapQ 20"
fi

echo "==========${0}===============" >> ${LOG}

exit $RVAL
