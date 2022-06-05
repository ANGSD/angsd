# beagle file reader test
# 220529 isinaltinkaya
if [[ ! $# -eq 1 ]]
then
    echo "Must supply an angsd binary; will exit"
	exit 1
else
    ANGSD=${1}
fi

LOG=${0}.log
rm -f ${LOG}

echo $ANGSD

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
REFDIR=${TDIR}/ref
RESDIR=${TDIR}/test_results
rm -rv ${RESDIR}
mkdir -pv ${RESDIR}






${ANGSD} -doMaf 4 -beagle ${TDIR}/test1.beagle -fai ${REFDIR}/ref1.fa.fai -out ${RESDIR}/test1
RVAL=0;

# 1_site
if [[ -n $($(cmp <(fz ${RESDIR}/test1.mafs.gz) <(f ${TDIR}/angsd_results/test1.mafs)) 2>> ${LOG}) ]];then
    echo "--------------"
    echo "Problem with beagleReader test1"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=1
fi



${ANGSD} -doMaf 4 -beagle ${TDIR}/test2.beagle -fai ${REFDIR}/ref2.fa.fai -out ${RESDIR}/test2
# chr1_site
if [[ -n $($(cmp <(fz ${RESDIR}/test2.mafs.gz) <(f ${TDIR}/angsd_results/test2.mafs)) 2>> ${LOG}) ]];then
    echo "--------------"
    echo "Problem with beagleReader test2"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=2
fi

${ANGSD} -doMaf 4 -beagle ${TDIR}/test3.beagle -fai ${REFDIR}/ref3.fa.fai -out ${RESDIR}/test3
# chr_1_site
if [[ -n $($(cmp <(fzd ${RESDIR}/test3.mafs.gz 2) <(f ${TDIR}/angsd_results/test1.mafs)) 2>> ${LOG}) ]];then
    echo "--------------"
    echo "Problem with beagleReader test3"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=3
fi


${ANGSD} -doMaf 4 -beagle ${TDIR}/test4.beagle -fai ${REFDIR}/ref4.fa.fai -out ${RESDIR}/test4
# chr_1_1_site
if [[ -n $($(cmp <(fzd ${RESDIR}/test4.mafs.gz 3) <(f ${TDIR}/angsd_results/test1.mafs)) 2>> ${LOG}) ]];then
    echo "--------------"
    echo "Problem with beagleReader test4"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=4
fi


${ANGSD} -doMaf 4 -beagle ${TDIR}/test5.beagle -fai ${REFDIR}/ref1.fa.fai -out ${RESDIR}/test5
# 1_site but \t separated
if [[ -n $($(cmp <(fz ${RESDIR}/test5.mafs.gz) <(f ${TDIR}/angsd_results/test1.mafs)) 2>> ${LOG}) ]];then
    echo "--------------"
    echo "Problem with beagleReader test5"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=1
fi


echo
echo "==========$0==============="
echo
exit $RVAL
