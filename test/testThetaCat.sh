#!/bin/bash
# thetaStat cat test. The two samples in bamlikes/ are copied to 3 contigs,
# thetas are made for all contigs and for each contig on its own (-r).
# Catting the per contig files must give the same index, print and do_stat
# output as the file with all contigs. Files with the same contig or a
# different nChr must give an error.
if [[ ! $# -eq 1 ]]
then
	echo "Must supply the folder with angsd and misc/; will exit"
	exit 1
else
	WDIR=${1}
fi

ANGSD=${WDIR}/angsd
REALSFS=${WDIR}/misc/realSFS
THETASTAT=${WDIR}/misc/thetaStat

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${WDIR}" >> ${LOG}

TDIR=thetacat
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with thetaStat cat: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

#the same reads on contig 1, 2 and 3
for f in s1 s2; do
	{
		printf '@HD\tVN:1.0\tSO:coordinate\n'
		for c in 1 2 3; do printf '@SQ\tSN:%s\tLN:150\n' ${c}; done
		for c in 1 2 3; do
			grep -v '^@' bamlikes/${f}.sam | awk -v c=${c} 'BEGIN{OFS="\t"}{$1=$1"_"c;$3=c;print}'
		done
	} > ${RESDIR}/${f}.sam
	echo ${RESDIR}/${f}.sam
done > ${RESDIR}/bam.filelist
head -1 ${RESDIR}/bam.filelist > ${RESDIR}/one.filelist

${ANGSD} -bam ${RESDIR}/bam.filelist -anc ${TDIR}/ref.fa -doSaf 1 -GL 1 -out ${RESDIR}/all >> ${LOG} 2>&1
${REALSFS} ${RESDIR}/all.saf.idx > ${RESDIR}/all.sfs 2>> ${LOG}
${REALSFS} saf2theta ${RESDIR}/all.saf.idx -sfs ${RESDIR}/all.sfs -outname ${RESDIR}/all >> ${LOG} 2>&1
for c in 1 2 3; do
	${REALSFS} saf2theta ${RESDIR}/all.saf.idx -sfs ${RESDIR}/all.sfs -r ${c} -outname ${RESDIR}/c${c} >> ${LOG} 2>&1
done
if [[ ! -s ${RESDIR}/all.thetas.idx || ! -s ${RESDIR}/c3.thetas.idx ]]; then
	problem "could not make the thetas files"
	exit 1
fi

#two files as arguments and one with -b
echo ${RESDIR}/c3.thetas.idx > ${RESDIR}/cat.filelist
if ! ${THETASTAT} cat ${RESDIR}/c1.thetas.idx ${RESDIR}/c2.thetas.idx -b ${RESDIR}/cat.filelist -outnames ${RESDIR}/merged >> ${LOG} 2>&1; then
	problem "cat exited with an error"
fi
if ! cmp ${RESDIR}/all.thetas.idx ${RESDIR}/merged.thetas.idx >> ${LOG} 2>&1; then
	problem "the index differs"
fi

for f in all merged; do
	${THETASTAT} print ${RESDIR}/${f}.thetas.idx > ${RESDIR}/${f}.print 2>> ${LOG}
	${THETASTAT} print ${RESDIR}/${f}.thetas.idx -r 2 > ${RESDIR}/${f}.print_r2 2>> ${LOG}
	${THETASTAT} do_stat ${RESDIR}/${f}.thetas.idx -outnames ${RESDIR}/${f}.whole >> ${LOG} 2>&1
	${THETASTAT} do_stat ${RESDIR}/${f}.thetas.idx -win 50 -step 10 -outnames ${RESDIR}/${f}.win >> ${LOG} 2>&1
	${THETASTAT} do_stat ${RESDIR}/${f}.thetas.idx -win 50 -step 10 -r 3 -outnames ${RESDIR}/${f}.win_r3 >> ${LOG} 2>&1
done
for x in print print_r2 whole.pestPG win.pestPG win_r3.pestPG; do
	if [[ ! -s ${RESDIR}/all.${x} ]] || ! cmp ${RESDIR}/all.${x} ${RESDIR}/merged.${x} >> ${LOG} 2>&1; then
		problem "${x} differs"
	fi
done

#the same contig twice
if ${THETASTAT} cat ${RESDIR}/c1.thetas.idx ${RESDIR}/all.thetas.idx -outnames ${RESDIR}/dup >> ${LOG} 2>&1; then
	problem "no error for the same contig in two files"
fi

#a different nChr (one sample)
${ANGSD} -bam ${RESDIR}/one.filelist -anc ${TDIR}/ref.fa -doSaf 1 -GL 1 -out ${RESDIR}/one >> ${LOG} 2>&1
${REALSFS} ${RESDIR}/one.saf.idx > ${RESDIR}/one.sfs 2>> ${LOG}
${REALSFS} saf2theta ${RESDIR}/one.saf.idx -sfs ${RESDIR}/one.sfs -r 3 -outname ${RESDIR}/one >> ${LOG} 2>&1
if [[ ! -s ${RESDIR}/one.thetas.idx ]]; then
	problem "could not make the thetas file with one sample"
fi
if ${THETASTAT} cat ${RESDIR}/c1.thetas.idx ${RESDIR}/one.thetas.idx -outnames ${RESDIR}/nchr >> ${LOG} 2>&1; then
	problem "no error for a different nChr"
fi

echo "==========$0===============" >> ${LOG}

exit $RVAL
