if [[ ! $# -eq 2 ]]
then
    echo "Supply an angsd binary and ASSODIR"
    exit 1
exit
fi

LOG=${0}.log
rm -f ${LOG}

ANGSD=$1
ASSODIR=$2

echo $ANGSD $ASSODIR


asso2new=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 2 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`

asso4new=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 4 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`

asso5new=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 5 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`

asso6new=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 6 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`


asso2binnew=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 2 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`

asso4binnew=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 4 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`

asso5binnew=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 5 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`

asso6binnew=`$ANGSD -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 6 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}|tail -n +1 | cut -f7`


## check results LRT and beta for each of the 8 analyses give back RVAL which one that fails

echo $asso2new

asso2=1.607704
RVAL=0
if [ ! "$asso2" = "$asso2new"  ] ;then
    echo "--------------"
    echo "Problem with score test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=2
fi

asso4=1.615975

if [ ! "$asso4" = "$asso4new"  ] ;then
    echo "--------------"
    echo "Problem with latent genotype test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=3
fi

asso5=1.607704

if [ ! "$asso5" = "$asso5new"  ] ;then
    echo "--------------"
    echo "Problem with hybrid test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=4
fi

asso6=1.615975

if [ ! "$asso6" = "$asso6new"  ] ;then
    echo "--------------"
    echo "Problem with dosage test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=5
fi

asso2bin=

if [ ! "$asso2bin" = "$asso2newbin"  ] ;then
    echo "--------------"
    echo "Problem with binary score test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=6
fi

asso4bin=

if [ ! "$asso4bin" = "$asso4newbin"  ] ;then
    echo "--------------"
    echo "Problem with binary latent genotype test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=7
fi

asso5bin=

if [ ! "$asso5bin" = "$asso5newbin"  ] ;then
    echo "--------------"
    echo "Problem with binary hybrid test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=8
fi

asso6bin=

if [ ! "$asso6bin" = "$asso6newbin"  ] ;then
    echo "--------------"
    echo "Problem with binary dosage test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=9
fi


exit $RVAL
