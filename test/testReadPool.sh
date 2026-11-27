#!/bin/bash
# The read pools in makeReadPool.cpp.
# When a file had no reads to give (an empty file, all reads filtered, or no
# reads in the -r region) and it was not the largest file, the end of that
# file was counted twice. With few files the chromosome/region was then
# stopped while the other files still had reads: most of the data was lost,
# or the rest of the chromosome came after the next one with some sites
# printed twice, all with exit code 0.
# Adding a file without reads must give the same sites and depths as the
# other files alone.
# readpool/A.bam has 100bp reads every 50bp on c and d, readpool/noc.bam a
# few reads on d only (made with samtools from the reads made below).
#
# ./testReadPool.sh angsd
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

TDIR=readpool
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with the read pools: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

hdr(){
	printf "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:30000\n@SQ\tSN:d\tLN:30000\n"
}
#reads chr from to step name mapq
reads(){
	awk -v chr=$1 -v from=$2 -v to=$3 -v step=$4 -v name=$5 -v mapq=$6 'BEGIN{OFS="\t"; s=sprintf("%100s",""); q=s; gsub(/ /,"A",s); gsub(/ /,"I",q); for(p=from;p<=to;p+=step) print name"_"chr"_"p,0,chr,p,mapq,"100M","*",0,0,s,q}'
}

{ hdr; reads c 1 20000 50 A 60; reads d 1 20000 50 A 60; } > ${RESDIR}/A.sam
hdr > ${RESDIR}/empty.sam
{ hdr; reads c 1 20000 2000 Z 0; } > ${RESDIR}/mapq0.sam

ARGS="-doCounts 1 -dumpCounts 1"
#name, files, extra options, the files of the expected output (same options), number of sites
TESTS=("empty_file"   "${RESDIR}/A.sam ${RESDIR}/empty.sam"   ""            "${RESDIR}/A.sam" 40100
       "all_filtered" "${RESDIR}/A.sam ${RESDIR}/mapq0.sam"   "-minMapQ 1"  "${RESDIR}/A.sam" 40100
       "no_reads_r"   "${TDIR}/A.bam ${TDIR}/noc.bam"         "-r c:"       "${TDIR}/A.bam"   20050
       "no_reads_rf"  "${TDIR}/A.bam ${TDIR}/noc.bam"         "-rf ${RESDIR}/rf.txt" "${TDIR}/A.bam" 20050)
printf "c:1-10000\nc:10001-30000\n" > ${RESDIR}/rf.txt

for ((i=0;i<${#TESTS[@]};i+=5)); do
	name=${TESTS[i]}
	for f in ${TESTS[i+1]}; do echo ${f}; done > ${RESDIR}/${name}.list
	for f in ${TESTS[i+3]}; do echo ${f}; done > ${RESDIR}/${name}.exp.list
	for t in ${name} ${name}.exp; do
		if ! ${ANGSD} -bam ${RESDIR}/${t}.list ${ARGS} ${TESTS[i+2]} -out ${RESDIR}/${t} >> ${LOG} 2>&1; then
			problem "${name}: '${TESTS[i+1]} ${TESTS[i+2]}' exited with an error"
		fi
	done
	gunzip -c ${RESDIR}/${name}.pos.gz > ${RESDIR}/${name}.pos
	gunzip -c ${RESDIR}/${name}.exp.pos.gz > ${RESDIR}/${name}.exp.pos
	if ! cmp -s ${RESDIR}/${name}.pos ${RESDIR}/${name}.exp.pos; then
		problem "${name}: the sites and depths of '${TESTS[i+1]} ${TESTS[i+2]}' are not the ones of '${TESTS[i+3]}' ($(tail -n +2 ${RESDIR}/${name}.pos | wc -l) vs $(tail -n +2 ${RESDIR}/${name}.exp.pos | wc -l) sites)"
	fi
	if [[ $(tail -n +2 ${RESDIR}/${name}.exp.pos | wc -l) -ne ${TESTS[i+4]} ]]; then
		problem "${name}: '${TESTS[i+3]}' alone does not have the ${TESTS[i+4]} sites of the reads"
	fi
done

#files with different headers must stop with exit code 1
{ printf "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:30000\n"; reads c 1 1000 50 C 60; } > ${RESDIR}/otherhdr.sam
printf "${RESDIR}/A.sam\n${RESDIR}/otherhdr.sam\n" > ${RESDIR}/otherhdr.list
${ANGSD} -bam ${RESDIR}/otherhdr.list ${ARGS} -out ${RESDIR}/otherhdr >> ${LOG} 2>&1
ex=$?
if [[ ${ex} -ne 1 ]]; then
	problem "files with different headers exited with ${ex}, not 1"
fi

echo "==========$0===============" >> ${LOG}

exit ${RVAL}
