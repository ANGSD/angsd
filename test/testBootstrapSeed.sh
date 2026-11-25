#!/bin/bash
# realSFS -bootstrap with the old default -seed -1 used the same seed in every run,
# so all runs gave the same bootstrap replicates, see #407. Now the default seed is
# from the clock. Two runs without -seed must differ, two runs with -seed 7 or
# -seed -1 (flat start values) must not, and -seed -1 must warn with -bootstrap.
#
# ./testBootstrapSeed.sh angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}
REALSFS=$(dirname ${ANGSD})/misc/realSFS

LOG=${0}.log
rm -f ${LOG}

RESDIR=bootstrapseed/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

${ANGSD} -bam bamlikes/bam.filelist -GL 1 -doSaf 1 -anc thetacat/ref.fa -out ${RESDIR}/pop >> ${LOG} 2>&1 || { echo "angsd -doSaf failed"; cat ${LOG}; exit 1; }

for r in a b; do
	${REALSFS} ${RESDIR}/pop.saf.idx -bootstrap 3 > ${RESDIR}/default_${r}.sfs 2>> ${LOG} || { echo "realSFS -bootstrap failed"; cat ${LOG}; exit 1; }
	sleep 1
	${REALSFS} ${RESDIR}/pop.saf.idx -bootstrap 3 -seed 7 > ${RESDIR}/seed7_${r}.sfs 2>> ${LOG} || { echo "realSFS -bootstrap -seed 7 failed"; cat ${LOG}; exit 1; }
	${REALSFS} ${RESDIR}/pop.saf.idx -bootstrap 3 -seed -1 > ${RESDIR}/flat_${r}.sfs 2> ${RESDIR}/flat_${r}.log || { echo "realSFS -bootstrap -seed -1 failed"; cat ${RESDIR}/flat_${r}.log; exit 1; }
	cat ${RESDIR}/flat_${r}.log >> ${LOG}
done

if [[ $(wc -l < ${RESDIR}/default_a.sfs) -ne 3 ]]; then
	echo "Problem with realSFS -bootstrap: expected 3 replicates in ${RESDIR}/default_a.sfs"
	exit 1
fi
if cmp -s ${RESDIR}/default_a.sfs ${RESDIR}/default_b.sfs; then
	echo "Problem with realSFS -bootstrap: two runs with the default seed gave the same replicates"
	exit 1
fi
if ! cmp -s ${RESDIR}/seed7_a.sfs ${RESDIR}/seed7_b.sfs; then
	echo "Problem with realSFS -bootstrap: two runs with -seed 7 gave different replicates"
	exit 1
fi
if ! cmp -s ${RESDIR}/flat_a.sfs ${RESDIR}/flat_b.sfs; then
	echo "Problem with realSFS -bootstrap: two runs with -seed -1 gave different replicates"
	exit 1
fi
if ! grep -q "same bootstrap replicates" ${RESDIR}/flat_a.log; then
	echo "Problem with realSFS -bootstrap: no warning with -seed -1"
	exit 1
fi
exit 0
