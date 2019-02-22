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


$ANGSD -vcf-gl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/en  >>${LOG} 2>&1
$ANGSD -vcf-gl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/to -r 1  >>${LOG} 2>&1
$ANGSD -vcf-gl ${VCF} -domajorminor 1 -domaf 1 -out ${ODIR}/tre -r 10  >>${LOG} 2>&1

md5sum  -c ${0}.md5>>${LOG} 2>&1 || exit 1
