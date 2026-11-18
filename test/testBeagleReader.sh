#!/bin/bash
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

function f(){
	cut -f1,3- $1|sed 1d
}

function fz(){
	gunzip -c $1 | cut -f1,3- |sed 1d
}


#handle delimiters
function fzd(){
	gunzip -c $1 | cut -f1,3- |sed 1d| cut -d_ -f${2}
}



echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD}" >> ${LOG}




TDIR=beagle_reader
REFDIR=${TDIR}/ref
RESDIR=${TDIR}/test_results
rm -rfv ${RESDIR} >> ${LOG}
mkdir -pv ${RESDIR} >> ${LOG}

${ANGSD} -doMaf 4 -beagle ${TDIR}/test1.beagle -fai ${REFDIR}/ref1.fa.fai -out ${RESDIR}/test1 >> ${LOG} 2>&1 
RVAL=0;

# 1_site
if ! cmp <(fz ${RESDIR}/test1.mafs.gz) <(f ${TDIR}/angsd_results/test1.mafs) >> ${LOG} 2>&1;then
    echo "--------------"
    echo "Problem with beagleReader test1"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=1
fi



${ANGSD} -doMaf 4 -beagle ${TDIR}/test2.beagle -fai ${REFDIR}/ref2.fa.fai -out ${RESDIR}/test2 >> ${LOG} 2>&1 
# chr1_site
if ! cmp <(fz ${RESDIR}/test2.mafs.gz) <(f ${TDIR}/angsd_results/test2.mafs) >> ${LOG} 2>&1;then
    echo "--------------"
    echo "Problem with beagleReader test2"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=2
fi

${ANGSD} -doMaf 4 -beagle ${TDIR}/test3.beagle -fai ${REFDIR}/ref3.fa.fai -out ${RESDIR}/test3 >> ${LOG} 2>&1 
# chr_1_site
if ! cmp <(fzd ${RESDIR}/test3.mafs.gz 2) <(f ${TDIR}/angsd_results/test1.mafs) >> ${LOG} 2>&1;then
    echo "--------------"
    echo "Problem with beagleReader test3"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=3
fi


${ANGSD} -doMaf 4 -beagle ${TDIR}/test4.beagle -fai ${REFDIR}/ref4.fa.fai -out ${RESDIR}/test4 >> ${LOG} 2>&1 
# chr_1_1_site
if ! cmp <(fzd ${RESDIR}/test4.mafs.gz 3) <(f ${TDIR}/angsd_results/test1.mafs) >> ${LOG} 2>&1;then
    echo "--------------"
    echo "Problem with beagleReader test4"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=4
fi


${ANGSD} -doMaf 4 -beagle ${TDIR}/test5.beagle -fai ${REFDIR}/ref1.fa.fai -out ${RESDIR}/test5 >> ${LOG} 2>&1 
# 1_site but \t separated
if ! cmp <(fz ${RESDIR}/test5.mafs.gz) <(f ${TDIR}/angsd_results/test1.mafs) >> ${LOG} 2>&1;then
    echo "--------------"
    echo "Problem with beagleReader test5"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=1
fi

${ANGSD} -doMaf 4 -beagle ${TDIR}/test6.beagle -fai ${REFDIR}/ref6.fa.fai -out ${RESDIR}/test6 >> ${LOG} 2>&1
# test1 data with long contig names and long numbers, so data lines are
# longer than the header; the reader used to overflow its line buffer on
# such lines (only reliably visible when built with -fsanitize=address)
if ! cmp <(fzd ${RESDIR}/test6.mafs.gz 11) <(f ${TDIR}/angsd_results/test1.mafs) >> ${LOG} 2>&1;then
    echo "--------------"
    echo "Problem with beagleReader test6"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=6
fi

echo "==========$0===============" >> ${LOG}

exit $RVAL
