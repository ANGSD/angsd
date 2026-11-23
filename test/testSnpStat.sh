#!/bin/bash
# -doSnpStat test.
# 1) 12 of the smallBam files (only when BAMDIR is given), the SNPs in
#    1:14000000-14030000 run with each filter. The .snpStat and .mafs output
#    must match snpstat/expected/, which was made with
#    angsd 0.943-10-g1c624c2 (before the snpstat fixes).
# 2) Without -doHWE the run must work and the HWE column is nan:nan.
# 3) Two samples with 1000bp reads made from snpstat/long_reads.txt, the reads
#    with the minor allele at c:1501 have MAPQ 255. The edge and mapQ Z scores
#    must be the ones from a Mann-Whitney test done outside angsd (distances
#    to the read end above 255 put in one bin).
# 4) Input without reads, options without what they need (the filters without
#    -doSnpStat 1, issue #401, and -SNP_pval with one individual, issue #649)
#    and -doSnpStat other than 0/1 must exit with 1.
#
# ./testSnpStat.sh angsd [BAMDIR]            runs the test
# ./testSnpStat.sh angsd BAMDIR expected     remakes snpstat/expected/ with this angsd
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}
BAMDIR=${2}
MAKEEXP=0
if [[ $# -eq 3 && "$3" == "expected" ]]; then
	MAKEEXP=1
fi

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD} and BAMDIR: '${BAMDIR}'" >> ${LOG}

TDIR=snpstat
EXPDIR=${TDIR}/expected
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with -doSnpStat: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

#1) the filters on smallBam
if [[ -n "${BAMDIR}" ]]; then
	ls ${BAMDIR}/*.bam | sort | head -12 > ${RESDIR}/bam.filelist
	BASE="-bam ${RESDIR}/bam.filelist -r 1:14000000-14030000 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doCounts 1 -doSnpStat 1"
	NAMES=(hwe sb edge hwefilter qscore mapq hetbias)
	ARGS=("-doHWE 1"
	      "-doHWE 1 -sb_pval 0.05"
	      "-doHWE 1 -edge_pval 0.5"
	      "-doHWE 1 -hwe_pval 0.05"
	      "-doHWE 1 -qscore_pval 0.05"
	      "-doHWE 1 -mapQ_pval 0.05"
	      "-doHWE 1 -doPost 1 -doGeno 1 -hetbias_pval 0.05")
	for i in ${!NAMES[@]}; do
		n=${NAMES[$i]}
		if ! ${ANGSD} ${BASE} ${ARGS[$i]} -out ${RESDIR}/${n} >> ${LOG} 2>&1; then
			problem "'${ARGS[$i]}' exited with an error"
			continue
		fi
		for t in snpStat mafs; do
			gunzip -c ${RESDIR}/${n}.${t}.gz > ${RESDIR}/${n}.${t}
			if [[ ${MAKEEXP} -eq 1 ]]; then
				cp ${RESDIR}/${n}.${t} ${EXPDIR}/${n}.${t}
			elif ! ./numcmp.sh ${EXPDIR}/${n}.${t} ${RESDIR}/${n}.${t} >> ${LOG} 2>&1; then
				problem "the .${t} for '${ARGS[$i]}' differs from ${EXPDIR}/${n}.${t}"
			fi
		done
	done
	if [[ ${MAKEEXP} -eq 1 ]]; then
		echo "Made the files in ${EXPDIR}"
		exit 0
	fi
else
	echo "No BAMDIR, the smallBam part is skipped" >> ${LOG}
fi

#2) without -doHWE the HWE column is nan:nan, with smallBam the other columns are as with -doHWE
if [[ -n "${BAMDIR}" ]]; then
	NOHWE="${BASE}"
else
	NOHWE="-bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doSnpStat 1"
fi
if ! ${ANGSD} ${NOHWE} -out ${RESDIR}/nohwe >> ${LOG} 2>&1; then
	problem "no -doHWE exited with an error"
else
	gunzip -c ${RESDIR}/nohwe.snpStat.gz > ${RESDIR}/nohwe.snpStat
	if [[ $(tail -n +2 ${RESDIR}/nohwe.snpStat | cut -f5 | sort -u) != "nan:nan" ]]; then
		problem "the HWE column is not nan:nan without -doHWE"
	fi
	if [[ -n "${BAMDIR}" ]]; then
		cut -f1-4,6- ${EXPDIR}/hwe.snpStat > ${RESDIR}/hwe.nohwecols
		cut -f1-4,6- ${RESDIR}/nohwe.snpStat > ${RESDIR}/nohwe.nohwecols
		if ! ./numcmp.sh ${RESDIR}/hwe.nohwecols ${RESDIR}/nohwe.nohwecols >> ${LOG} 2>&1; then
			problem "the output without -doHWE differs in other columns than HWE"
		fi
	fi
fi

#3) long reads and MAPQ 255. The reads are the reference, the ones with alt=1 have the next base at c:1501
for s in s1 s2; do
	awk -v s=${s} -v snp=1501 'BEGIN{OFS="\t"}
		NR==FNR{if($0!~/^>/)ref=ref $0; next}
		/^#/{next}
		index($1,s"_")==1{
			seq=substr(ref,$3,1000)
			if($5==1){
				b=substr(ref,snp,1); a=(b=="A")?"C":(b=="C")?"G":(b=="G")?"T":"A"
				o=snp-$3+1; seq=substr(seq,1,o-1) a substr(seq,o+1)
			}
			q=sprintf("%1000s",""); gsub(/ /,"I",q)
			print $1,$2,"c",$3,$4,"1000M","*",0,0,seq,q
		}' ${TDIR}/long_ref.fa ${TDIR}/long_reads.txt | { printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:c\tLN:3000\n'; cat; } > ${RESDIR}/long_${s}.sam
	echo ${RESDIR}/long_${s}.sam
done > ${RESDIR}/long.filelist
if ! ${ANGSD} -bam ${RESDIR}/long.filelist -ref ${TDIR}/long_ref.fa -GL 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doHWE 1 -doSnpStat 1 -out ${RESDIR}/long >> ${LOG} 2>&1; then
	problem "the long reads exited with an error"
else
	line=$(gunzip -c ${RESDIR}/long.snpStat.gz | awk -F'\t' '$2==1501')
	mapqZ=$(echo "${line}" | cut -f7 | cut -d: -f1)
	edgeZ=$(echo "${line}" | cut -f8 | cut -d: -f1)
	echo "long reads c:1501 mapQ Z: ${mapqZ} edge Z: ${edgeZ}" >> ${LOG}
	if ! awk -v a="${mapqZ}" 'BEGIN{d=a-(-7.698004);exit !(a!="" && d<1e-5 && d>-1e-5)}'; then
		problem "mapQ Z with MAPQ 255 is '${mapqZ}', should be -7.698004"
	fi
	if ! awk -v a="${edgeZ}" 'BEGIN{d=a-(-0.923760);exit !(a!="" && d<1e-5 && d>-1e-5)}'; then
		problem "edge Z with 1000bp reads is '${edgeZ}', should be -0.923760"
	fi
fi

#4) errors
GLF=tajima/output/glout.glf.gz
if [[ ! -f ${GLF} ]]; then
	GLF=${RESDIR}/empty.glf.gz
	printf '' | gzip > ${GLF}
fi
SMALL="-bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doSnpStat 1"
#-SNP_pval needs more than one individual (#649)
echo bamlikes/s1.sam > ${RESDIR}/one.filelist
ERRS=("-bam ${RESDIR}/one.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6"
      "-glf ${GLF} -nInd 20 -fai hg19.fa.fai -isSim 1 -doMajorMinor 1 -doMaf 1 -doHWE 1 -doSnpStat 1"
      "${SMALL} -hwe_pval 0.05"
      "${SMALL} -doHWE 1 -hetbias_pval 0.05"
      "-bam bamlikes/bam.filelist -GL 1 -doSnpStat 1"
      "-bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -doSnpStat 2"
      "-bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -sb_pval 0.05"
      "-bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -doSnpStat 0 -qscore_pval 0.05")
for e in "${ERRS[@]}"; do
	${ANGSD} ${e} -out ${RESDIR}/err >> ${LOG} 2>&1
	ex=$?
	if [[ ${ex} -ne 1 ]]; then
		problem "exit code ${ex}, should be 1, for: ${e}"
	fi
done

#the filters without -doSnpStat 1 used to be reported as unknown arguments
msg=$(${ANGSD} -bam bamlikes/bam.filelist -GL 1 -doMajorMinor 1 -doMaf 1 -sb_pval 0.05 -out ${RESDIR}/err 2>&1)
if ! grep -q -- "-sb_pval needs -doSnpStat 1" <<< "${msg}"; then
	problem "-sb_pval without -doSnpStat 1 does not say it needs -doSnpStat 1"
fi

echo "==========$0===============" >> ${LOG}

exit $RVAL
