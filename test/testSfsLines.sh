#!/bin/bash
# realSFS -nSites writes one sfs per block of sites, one per line. Such a file
# given with -sfs used to stop with 'Pxroblem with size of dimension of prior'
# and exit code 0 (#530). It must now say that the file has several lines,
# exit with 1, and work when the lines are summed. A file with the wrong number
# of values or no values must also exit with 1, and so must saf files that
# can't be opened (#646).
# Runs on the two samples in bamlikes/ with -doSaf 1 (thetacat/ref.fa as -anc).
#
# ./testSfsLines.sh angsd
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
	echo "Problem with realSFS -sfs files: no ${REALSFS}"
	exit 1
fi

RESDIR=sfslines/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with realSFS -sfs files: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

if ! ${ANGSD} -bam bamlikes/bam.filelist -GL 1 -doSaf 1 -anc thetacat/ref.fa -out ${RESDIR}/s >> ${LOG} 2>&1; then
	problem "-doSaf exited with an error"
	exit 1
fi
SAF=${RESDIR}/s.saf.idx

#blocks of 40 sites, one sfs per line
if ! ${REALSFS} ${SAF} -nSites 40 -maxIter 5 > ${RESDIR}/blocks.sfs 2>>${LOG}; then
	problem "realSFS -nSites exited with an error"
	exit 1
fi
if [[ $(wc -l < ${RESDIR}/blocks.sfs) -lt 2 ]]; then
	problem "realSFS -nSites 40 gave less than 2 lines"
	exit 1
fi
awk '{for(i=1;i<=NF;i++)s[i]+=$i}END{for(i=1;i<=NF;i++)printf "%f ",s[i];print ""}' ${RESDIR}/blocks.sfs > ${RESDIR}/sum.sfs
awk '{print $1,$2,$3}' ${RESDIR}/sum.sfs > ${RESDIR}/short.sfs
: > ${RESDIR}/empty.sfs

#saf2theta and the sfs optimisation both read -sfs
for prog in saf2theta optim; do
	if [[ ${prog} == saf2theta ]]; then
		run(){ ${REALSFS} saf2theta ${SAF} -sfs $1 -outname ${RESDIR}/theta_$2; }
	else
		run(){ ${REALSFS} ${SAF} -sfs $1 -maxIter 5 > ${RESDIR}/optim_$2.sfs; }
	fi

	run ${RESDIR}/blocks.sfs blocks > ${RESDIR}/${prog}_blocks.log 2>&1
	e=$?
	cat ${RESDIR}/${prog}_blocks.log >> ${LOG}
	if [[ ${e} -ne 1 ]]; then
		problem "${prog} with an sfs on several lines exited with ${e}, not 1"
	fi
	if ! grep -q "lines with" ${RESDIR}/${prog}_blocks.log; then
		problem "${prog} with an sfs on several lines doesn't say the file has several lines"
	fi

	if ! run ${RESDIR}/sum.sfs sum >> ${LOG} 2>&1; then
		problem "${prog} with the summed sfs exited with an error"
	fi

	for bad in short empty; do
		run ${RESDIR}/${bad}.sfs ${bad} >> ${LOG} 2>&1
		e=$?
		if [[ ${e} -ne 1 ]]; then
			problem "${prog} with a ${bad} sfs file exited with ${e}, not 1"
		fi
	done
done

#a saf file that can't be opened must exit with 1, it used to exit with 0 (#646)
cp ${SAF} ${RESDIR}/nogz.saf.idx
for idx in ${RESDIR}/missing.saf.idx ${RESDIR}/nogz.saf.idx; do
	${REALSFS} ${idx} > /dev/null 2>>${LOG}
	e=$?
	if [[ ${e} -ne 1 ]]; then
		problem "realSFS with ${idx} (missing saf files) exited with ${e}, not 1"
	fi
	${REALSFS} saf2theta ${idx} -sfs ${RESDIR}/sum.sfs -outname ${RESDIR}/theta_missing > /dev/null 2>>${LOG}
	e=$?
	if [[ ${e} -ne 1 ]]; then
		problem "saf2theta with ${idx} (missing saf files) exited with ${e}, not 1"
	fi
done

exit ${RVAL}
