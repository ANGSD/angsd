# beagle file reader test
# 220529 isinaltinkaya
if [[ ! $# -eq 1 ]]
then
    echo "Must supply an angsd binary; will exit"
	exit 1
else
    ANGSD=${1}
fi


function f(){
	cut -f2 --complement $1|sed 1d
}

function fz(){
	gunzip -c $1 | cut -f2 --complement |sed 1d
}


#handle delimiters
function fzd(){
	gunzip -c $1 | cut -f2 --complement |sed 1d| cut -d_ -f${2}
}



echo "==========${0}=========="
echo "Running test: ${0} with ${ANGSD}"
echo



TDIR=beagle_reader
REFDIR=${TDIR}/ref/
RESDIR=${TDIR}/test_results/
rm -rv ${RESDIR}
mkdir -pv ${RESDIR}

${ANGSD} -doMaf 4 -beagle ${TDIR}/test1.beagle -fai ${REFDIR}/ref1.fa.fai -out ${RESDIR}/test1
${ANGSD} -doMaf 4 -beagle ${TDIR}/test2.beagle -fai ${REFDIR}/ref2.fa.fai -out ${RESDIR}/test2
${ANGSD} -doMaf 4 -beagle ${TDIR}/test3.beagle -fai ${REFDIR}/ref3.fa.fai -out ${RESDIR}/test3
${ANGSD} -doMaf 4 -beagle ${TDIR}/test4.beagle -fai ${REFDIR}/ref4.fa.fai -out ${RESDIR}/test4

# 1_site
cmp <(fz ${RESDIR}/test1.mafs.gz) <(f ${TDIR}/angsd_results/test1.mafs)

# chr1_site
cmp <(fz ${RESDIR}/test2.mafs.gz) <(f ${TDIR}/angsd_results/test2.mafs)

# chr_1_site
cmp <(fzd ${RESDIR}/test3.mafs.gz 2) <(f ${TDIR}/angsd_results/test1.mafs)

# chr_1_1_site
cmp <(fzd ${RESDIR}/test4.mafs.gz 3) <(f ${TDIR}/angsd_results/test1.mafs)

echo
echo "==========$0==============="
echo
