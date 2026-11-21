if [[ ! $# -eq 2 ]]
then
    echo "Supply an angsd binary and BAMdir"
    exit 1
fi

LOG=${0}.log

MD5=md5sum
if [[ "$OSTYPE" == "darwin"* ]]; then
    MD5="$(cd "$(dirname "$0")" && pwd)/md5osx.sh"
fi
rm -f ${LOG}

ANGSD=$1
BDIR=$2

echo $ANGSD $BDIR

ls $BDIR/*.bam >smallBam.filelist

#everything made by the test goes to sites/test_results, the sites files are copied
#there as 'angsd sites index' writes the .bin and .idx files next to them
RES=sites/test_results
rm -rf ${RES}
mkdir -p ${RES}
cp sites/s1 sites/s2 sites/s3 sites/s4 sites/s5 ${RES}/

#runs angsd and stops the test with exit code 7 if angsd fails
runangsd(){
    ${ANGSD} "$@" 2>>${LOG}
    ret=$?
    if [[ ${ret} -ne 0 ]]; then
	echo "Problem: angsd $* gave exit code ${ret}"
	exit 7
    fi
}

runangsd sites index ${RES}/s1
runangsd -out ${RES}/run1 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 1 -sites ${RES}/s1

runangsd sites index ${RES}/s2 -compl 1
runangsd -out ${RES}/run2 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 1 -sites ${RES}/s2

runangsd sites index ${RES}/s3
runangsd -out ${RES}/run3 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 1 -sites ${RES}/s3

runangsd sites index ${RES}/s4
runangsd -out ${RES}/run4 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 3 -sites ${RES}/s4

runangsd sites index ${RES}/s5
runangsd -out ${RES}/run5 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 3 -sites ${RES}/s5

runangsd -out ${RES}/run6 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 3 -sites ${RES}/s5 -dobcf 1 -dopost 1


cd ${RES}
${MD5}  -c ../md5orig.orig >>../../${LOG} 2>&1||exit 8

