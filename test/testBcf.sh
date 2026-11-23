#!/bin/bash
# -doBcf test of FORMAT/AD. With -doCounts the AD of each sample must be the
# number of bases of the major and the minor allele (REF,ALT) as counted by
# -dumpCounts 4 (.counts.gz) with the major and minor from .mafs.gz.
# Runs on the two samples in bamlikes/ and, when BAMDIR is given, on the SNPs
# of 12 of the smallBam files in 1:14000000-14030000.
# The bcf is read with htsfile from htslib.
#
# ./testBcf.sh angsd [BAMDIR]
if [[ $# -lt 1 ]]
then
	echo "Must supply an angsd binary; will exit"
	exit 1
fi
ANGSD=${1}
BAMDIR=${2}

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD} and BAMDIR: '${BAMDIR}'" >> ${LOG}

HTSFILE=$(dirname ${ANGSD})/htslib/htsfile
if [[ ! -x ${HTSFILE} ]]; then
	HTSFILE=$(command -v htsfile)
fi
if [[ -z "${HTSFILE}" ]]; then
	echo "No htsfile to read the bcf with, -doBcf is not tested" >> ${LOG}
	exit 0
fi

RESDIR=bcf/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with -doBcf: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

#checks AD in NAME.bcf against NAME.counts.gz and NAME.pos.gz. The alleles are
#major,minor from NAME.mafs.gz, or with MODE mcall the REF and ALT alleles of the record
checkad(){
	local n=$1
	local mode=$2
	${HTSFILE} -c ${RESDIR}/${n}.bcf > ${RESDIR}/${n}.vcf 2>> ${LOG}
	if ! grep -q '^##FORMAT=<ID=AD,' ${RESDIR}/${n}.vcf; then
		problem "${n}: no FORMAT/AD in the header"
	fi
	if grep -q 'ID=DPR,' ${RESDIR}/${n}.vcf; then
		problem "${n}: DPR is in the header but never written"
	fi
	paste <(gunzip -c ${RESDIR}/${n}.pos.gz) <(gunzip -c ${RESDIR}/${n}.counts.gz) | tail -n +2 > ${RESDIR}/${n}.poscounts
	gunzip -c ${RESDIR}/${n}.mafs.gz | tail -n +2 > ${RESDIR}/${n}.mafs
	if ! awk -F'\t' -v out=${RESDIR}/${n}.adcheck -v mode="${mode}" '
		BEGIN{b["A"]=0;b["C"]=1;b["G"]=2;b["T"]=3}
		FILENAME~/poscounts$/{k=$1":"$2; for(i=4;i<=NF;i++) cnt[k,i-4]=$i; next}
		FILENAME~/mafs$/{maj[$1":"$2]=$3; mi[$1":"$2]=$4; next}
		/^#/{next}
		{
			k=$1":"$2
			if(mode=="mcall"){
				nal=1; al[1]=$4
				if($5!=".") nal+=split($5,alt,",")
				for(i=2;i<=nal;i++) al[i]=alt[i-1]
			}else{
				if(!(k in maj)){print "no mafs line for "k; bad++; next}
				if($4!=maj[k]||$5!=mi[k]){print k": REF,ALT is "$4","$5", mafs has "maj[k]","mi[k]; bad++}
				nal=2; al[1]=maj[k]; al[2]=mi[k]
			}
			nf=split($9,f,":"); ad=0
			for(i=1;i<=nf;i++) if(f[i]=="AD") ad=i
			if(ad==0){print k": no AD in FORMAT "$9; bad++; next}
			if(nal>2) multi++
			for(j=10;j<=NF;j++){
				split($j,v,":"); s=j-10
				want=cnt[k,4*s+b[al[1]]]
				for(i=2;i<=nal;i++) want=want","cnt[k,4*s+b[al[i]]]
				if(v[ad]!=want){print k" sample "s": AD is "v[ad]", counts give "want; bad++}
				checked++
			}
		}
		END{print checked" sample AD values checked, "multi+0" records with more than one ALT" > out; if(checked==0) print "no records"; exit (bad>0||checked==0)}' ${RESDIR}/${n}.poscounts ${RESDIR}/${n}.mafs ${RESDIR}/${n}.vcf >> ${LOG} 2>&1; then
		problem "${n}: AD differs from the counts"
	fi
	cat ${RESDIR}/${n}.adcheck >> ${LOG} 2>/dev/null
}

ARGS="-GL 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 -dumpCounts 4 -doPost 1 -doGeno 1 -doBcf 1"

if ${ANGSD} -bam bamlikes/bam.filelist ${ARGS} -out ${RESDIR}/bamlikes >> ${LOG} 2>&1; then
	checkad bamlikes
	#a sample without reads has no data, GL, PL and GP must be missing and not 0,0,0 (issue 553)
	if ! awk -F'\t' '/^#/{next}
		{
			nf=split($9,f,":")
			for(i=1;i<=nf;i++) c[f[i]]=i
			for(j=10;j<=NF;j++){
				split($j,v,":")
				miss=(v[c["GL"]]==".,.,."&&v[c["PL"]]==".,.,."&&v[c["GP"]]==".,.,.")
				if((v[c["DP"]]==0)!=miss){print $1":"$2" sample "j-10": DP is "v[c["DP"]]" and GL:PL:GP is "v[c["GL"]]":"v[c["PL"]]":"v[c["GP"]]; bad++}
				if(v[c["DP"]]==0) nodata++
			}
		}
		END{if(nodata==0) print "no samples without reads"; exit (bad>0||nodata==0)}' ${RESDIR}/bamlikes.vcf >> ${LOG} 2>&1; then
		problem "bamlikes: GL/PL/GP of samples without reads are not missing"
	fi
else
	problem "bamlikes exited with an error"
fi

#-doMcall can give more than one ALT. It needs -ref, the random reference of the
#same length in thetacat/ is used, so many sites have reads with 2 other bases
if ${ANGSD} -bam bamlikes/bam.filelist -ref thetacat/ref.fa -doMcall 1 ${ARGS} -out ${RESDIR}/mcall >> ${LOG} 2>&1; then
	checkad mcall mcall
	if ! grep -q 'more than one ALT' ${RESDIR}/mcall.adcheck || grep -q ' 0 records with more than one ALT' ${RESDIR}/mcall.adcheck; then
		problem "-doMcall gave no records with more than one ALT, so they are not tested"
	fi
else
	problem "-doMcall exited with an error"
fi

if [[ -n "${BAMDIR}" ]]; then
	ls ${BAMDIR}/*.bam | sort | head -12 > ${RESDIR}/bam.filelist
	if ${ANGSD} -bam ${RESDIR}/bam.filelist -r 1:14000000-14030000 -SNP_pval 1e-6 ${ARGS} -out ${RESDIR}/smallbam >> ${LOG} 2>&1; then
		checkad smallbam
	else
		problem "smallBam exited with an error"
	fi
fi

echo "==========$0===============" >> ${LOG}

exit $RVAL
