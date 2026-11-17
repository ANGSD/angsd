if [[ ! $# -eq 1 ]]
then
    echo "Supply an angsd binary"
    exit 1
exit
fi

LOG=${0}.log

MD5=md5sum
if [[ "$OSTYPE" == "darwin"* ]]; then
    MD5="$(cd "$(dirname "$0")" && pwd)/md5osx.sh"
fi
rm -f ${LOG}

ANGSD=$1
echo $ANGSD $VCF

ODIR=${0}.dir
#rm -rf ${ODIR}
mkdir -p ${ODIR}


$ANGSD -b smallBam.filelist -docounts 1 -dohaplocall 1 -out ${ODIR}/haplo -nind 2  >>${LOG} 2>&1

${MD5}  -c ${0}.md5>>${LOG} 2>&1 || exit 1
