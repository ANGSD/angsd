if [[ ! $# -eq 1 ]]
then
    echo "Supply an angsd binary"
    exit 1
exit
fi

LOG=${0}.log
rm -f ${LOG}

ANGSD=$1
echo $ANGSD $VCF

ODIR=${0}.dir
#rm -rf ${ODIR}
mkdir -p ${ODIR}


$ANGSD -b smallBam.filelist -docounts 1 -dohaplocall 1 -out ${ODIR}/haplo -nind 2  >>${LOG} 2>&1

md5sum  -c ${0}.md5>>${LOG} 2>&1 || exit 1
