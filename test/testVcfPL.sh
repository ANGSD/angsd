#!/bin/bash
# vcf input with large PL values (issue 585). A PL above ~3080 was read as a
# likelihood of 0 (-inf), so with -SNP_pval the LRT of a site where a sample
# is homozygous for the minor allele with such a PL for 0/0 was NaN, and angsd
# exited with 0 and empty output. Also checks that -vcf-gl on a vcf with only
# PL exits with an error instead of skipping all sites.
#
# ./testVcfPL.sh angsd
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

RESDIR=vcfpl/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

RVAL=0

problem(){
	echo "--------------"
	echo "Problem with vcf PL input: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

VCF=${RESDIR}/in.vcf
printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">\n' > ${VCF}
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n' >> ${VCF}
printf '1\t10\t.\tG\tA\t.\tPASS\t.\tPL\t0,60,900\t0,90,1350\t3333,187,0\t0,0,0\n' >> ${VCF}
printf '1\t20\t.\tC\tT\t.\tPASS\t.\tPL\t0,60,900\t0,30,450\t0,90,1350\t0,45,600\n' >> ${VCF}

if ! ${ANGSD} -vcf-pl ${VCF} -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -out ${RESDIR}/pl >> ${LOG} 2>&1; then
	problem "-vcf-pl with -SNP_pval exited with an error"
fi
if ! gunzip -c ${RESDIR}/pl.mafs.gz 2>> ${LOG} | awk -F'\t' '$1=="1"&&$2==10{f=1}END{exit !f}'; then
	problem "the site with PL 3333 is not in the .mafs"
fi

if ${ANGSD} -vcf-gl ${VCF} -doMajorMinor 1 -doMaf 1 -out ${RESDIR}/gl >> ${LOG} 2>&1; then
	problem "-vcf-gl on a vcf without GL did not exit with an error"
fi

echo "==========${0}===============" >> ${LOG}

exit $RVAL
