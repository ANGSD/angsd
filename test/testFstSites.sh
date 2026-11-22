#!/bin/bash
# realSFS with -sites on saf files where some chromosomes have no -sites.
# fst index used another chromosome for those (written under the wrong name),
# and could stop with 'SAF file is out of sync' (exit code 0), see #627.
# realSFS -r chr -sites used the other chromosomes when chr had no -sites.
# Two pops (bamlikes/s1.sam and s2.sam) with the same reads on chr 1, 2 and 3
# of thetacat/ref.fa. fst index must give the sites in -sites, with the same
# values as without -sites, and realSFS must not read other chromosomes.
#
# ./testFstSites.sh angsd
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
	echo "Problem with realSFS and -sites: no ${REALSFS}"
	exit 1
fi

RESDIR=fstsites/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with realSFS and -sites: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

#the reads of one sample on chr 1, 2 and 3
for s in s1 s2; do
	{
		printf "@HD\tVN:1.0\tSO:coordinate\n"
		for c in 1 2 3; do printf "@SQ\tSN:${c}\tLN:150\n"; done
		for c in 1 2 3; do
			grep -v '^@' bamlikes/${s}.sam | awk -v c=${c} 'BEGIN{OFS="\t"}{$3=c;print}'
		done
	} > ${RESDIR}/${s}.sam
	echo ${RESDIR}/${s}.sam > ${RESDIR}/${s}.filelist
	if ! ${ANGSD} -bam ${RESDIR}/${s}.filelist -GL 1 -doSaf 1 -anc thetacat/ref.fa -out ${RESDIR}/${s} >> ${LOG} 2>&1; then
		problem "-doSaf on ${s} exited with an error"
		exit 1
	fi
done
P1=${RESDIR}/s1.saf.idx
P2=${RESDIR}/s2.saf.idx

${REALSFS} print ${P1} ${P2} 2>>${LOG} | cut -f1,2 > ${RESDIR}/shared.pos
for c in 1 2 3; do
	if [[ $(awk -v c=${c} '$1==c' ${RESDIR}/shared.pos | wc -l) -lt 10 ]]; then
		problem "less than 10 shared sites on chr ${c}"
		exit 1
	fi
done

if ! ${REALSFS} ${P1} ${P2} -maxIter 5 > ${RESDIR}/2d.sfs 2>>${LOG}; then
	problem "realSFS 2d sfs exited with an error"
fi
if ! ${REALSFS} fst index ${P1} ${P2} -sfs ${RESDIR}/2d.sfs -fstout ${RESDIR}/all >> ${LOG} 2>&1; then
	problem "fst index without -sites exited with an error"
fi
${REALSFS} fst print ${RESDIR}/all.fst.idx 2>>${LOG} | sort > ${RESDIR}/all.fst.txt

#name: which chromosomes have -sites
for name in 2 12 23 13; do
	sites=${RESDIR}/sites${name}
	awk -v n=${name} 'index(n,$1)&&NR%2' ${RESDIR}/shared.pos > ${sites}
	if ! ${ANGSD} sites index ${sites} >> ${LOG} 2>&1; then
		problem "sites index on ${sites} exited with an error"
		continue
	fi
	out=${RESDIR}/fst${name}
	if ! ${REALSFS} fst index ${P1} ${P2} -sfs ${RESDIR}/2d.sfs -fstout ${out} -sites ${sites} > ${out}.log 2>&1; then
		cat ${out}.log >> ${LOG}
		problem "fst index with -sites on chr ${name} exited with an error"
		continue
	fi
	cat ${out}.log >> ${LOG}
	if grep -q "out of sync" ${out}.log; then
		problem "fst index with -sites on chr ${name}: saf files out of sync"
	fi
	${REALSFS} fst print ${out}.fst.idx 2>>${LOG} | sort > ${out}.fst.txt
	if ! cmp -s <(cut -f1,2 ${out}.fst.txt) <(sort ${sites}); then
		problem "fst index with -sites on chr ${name} doesn't give the sites in -sites"
	fi
	#the values of a site don't depend on -sites
	if [[ -n $(join -t $'\t' -v1 <(awk '{print $1":"$2"\t"$3"\t"$4}' ${out}.fst.txt) <(awk '{print $1":"$2"\t"$3"\t"$4}' ${RESDIR}/all.fst.txt)) ]]; then
		problem "fst index with -sites on chr ${name} gives other values than without -sites"
	fi

	#2d sfs and print with -sites read all pops in sync
	if ! ${REALSFS} ${P1} ${P2} -sites ${sites} -maxIter 5 > ${out}.sfs 2> ${out}.sfs.log; then
		cat ${out}.sfs.log >> ${LOG}
		problem "2d sfs with -sites on chr ${name} exited with an error"
	fi
	if ! ${REALSFS} print ${P1} ${P2} -sites ${sites} 2>>${LOG} | cut -f1,2 | cmp -s - ${sites}; then
		problem "realSFS print with -sites on chr ${name} doesn't give the sites in -sites"
	fi
done

#-r on a chromosome without -sites must not read the other chromosomes
${REALSFS} print ${P1} ${P2} -r 3 -sites ${RESDIR}/sites2 > ${RESDIR}/r3.txt 2>>${LOG}
if [[ -s ${RESDIR}/r3.txt ]]; then
	problem "realSFS print -r 3 with -sites only on chr 2 gives sites"
fi
${REALSFS} ${P1} ${P2} -r 3 -sites ${RESDIR}/sites2 -maxIter 5 > ${RESDIR}/r3.sfs 2>>${LOG}
if [[ -s ${RESDIR}/r3.sfs ]]; then
	problem "realSFS -r 3 with -sites only on chr 2 gives an sfs"
fi

exit ${RVAL}
