#!/bin/bash
#run locally like
#make test BAMDIR=smallBam/

PRG=""
BAMDIR=""

if [ $# -eq 0 ] 
then
    exit 1;
fi

if [ $# -eq 1 ]
then
    PRG=$1
fi

if [ $# -eq 2 ]
then
    PRG=$1
    BAMDIR=$2
fi


echo "--------------------"
echo "Using PRG: '${PRG}' and BAMDIR: '${BAMDIR}'"
echo "--------------------"

WDIR=`dirname $PRG`

RVAL=0

if [[ ! -z "$BAMDIR" ]]; then
echo "Testing -sites"
time ./testFilterSites.sh $WDIR/angsd $BAMDIR
if [ ! $? -eq 0  ]   ;then
    echo "Problem with -sites exit code: $?"
    cat ./testFilterSites.sh.log
    RVAL=1
fi
fi

if [[ ! -z "$BAMDIR" ]]; then
echo "Testing vcfreading"
time ./testVcf.sh $WDIR/angsd ${BAMDIR}/small2.bcf
if [ ! $? -eq 0  ]   ;then
    echo "Problem with -vcf-pl exit code: $?"
    cat ./testVcf.sh.log
    RVAL=1
fi
fi

#echo "Testing neutrality test statistics"
time ./testTaj.sh $WDIR
if [ ! $? -eq 0 ] ;then
    echo "Problem with neutrality test statistics exit code: $?"
    cat ./testTaj.sh.log
    RVAL=1
fi

echo "Testing fst using msms"
time ./testFst.sh $WDIR
if [ ! $? -eq 0 ] ;then
    echo "Problem with Fst test statistics exit code: $?"
    cat ./testFst.sh.log
    RVAL=1
fi

echo "Testing fst_folded using msms"
time ./testFst_folded.sh $WDIR
if [ ! $? -eq 0 ] ;then
    echo "Problem with Fst_folded test statistics exit code: $?"
    cat ./testFst_folded.sh.log
    RVAL=1
fi

echo "Testing saf_beagle"
time ./testSaf_beagle.sh $WDIR
if [ ! $? -eq 0 ] ;then
    echo "Problem with saf_beagle test statistics exit code: $?"
    cat ./testSaf_beagle.sh.log
    RVAL=1
fi

echo "Testing SFS"
time ./testSFS.sh $WDIR
if [ ! $? -eq 0  ]   ;then
    echo "Problem with SFS exit code: $?"
    cat ./testSFS.sh.log
    RVAL=1
fi

if [[ ! -z "$BAMDIR" ]]; then
echo "Testing basic mpileup"
time ./testBam.sh $WDIR/angsd $BAMDIR
if [ ! $? -eq 0  ]   ;then
    echo "Problem with basic pileup exit code: $?"
    cat ./testBam.sh.log
    RVAL=1
fi
fi

echo "Testing association"
time ./testDoAsso2456.sh $WDIR
if [ ! $? -eq 0  ]   ;then
    echo "Problem with association exit code: $?"
    cat ./testDoAsso2456.sh.log
    RVAL=1
fi

if [[ ! -z "$BAMDIR" ]]; then
echo "Testing haplocall"
time ./testHaploCall.sh $WDIR/angsd
if [ ! $? -eq 0  ]   ;then
    echo "Problem with haplocall exit code: $?"
    cat ./testHaploCall.sh.log
    RVAL=1
fi
fi

echo "Testing beagle reader"
time ./testBeagleReader.sh $WDIR/angsd
if [ ! $? -eq 0  ]   ;then
    echo "Problem with beagle reader exit code: $?"
    RVAL=1
fi


exit ${RVAL}

if false; then
    
    echo testAsso6.sh $PRG
    ./testAsso6.sh $PRG
    
    #echo ./testErr.sh $PRG
    #./testErr.sh $PRG
    
    echo ./testGL6.sh $PRG
    ./testGL6.sh $PRG
    
    echo ./testMisc9.sh $PRG
    ./testMisc9.sh $PRG
    
    echo ./testAbba.sh $PRG
    ./testAbba.sh $PRG
        
    echo ./testFasta.sh
    ./testFasta.sh $PRG
    
    echo ./testBaq.sh $PRG
    ./testBaq3.sh $PRG
    
    echo "Netaccess is now deprecated"
    #echo ./testNetAccess.sh $PRG
    #./testNetAccess.sh $PRG
fi
