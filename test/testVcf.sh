if [[ ! $# -eq 2 ]]
then
    echo "Supply an angsd binary and bcf"
    exit 1
exit
fi

LOG=${0}.log
rm -f ${LOG}

ANGSD=$1
VCF=$2
echo $ANGSD $VCF

ODIR=${0}.dir
#rm -rf ${ODIR}
mkdir -p ${ODIR}


$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/en  >>${LOG} 2>&1
$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/to -r 1  >>${LOG} 2>&1
$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/tre -r 10  >>${LOG} 2>&1
$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/fire -r 1:14000303  >>${LOG} 2>&1
$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/fem -r 1:14000302  >>${LOG} 2>&1
$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/seks -r 7:14096608-14098164  >>${LOG} 2>&1
$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/syv -r 7:14096608-14098165  >>${LOG} 2>&1

$ANGSD -vcf-pl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/otte -rf sites.rf >>${LOG} 2>&1

md5sum  -c ${0}.md5 >>${LOG} 2>&1 || exit 1
