#!/bin/bash

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

echo "Testing neutrality test statistics"
#./testTaj.sh $WDIR
if [ ! $? -eq 0 ] ;then
    echo "Problem with neutrality test statistics exit code: $?"
    cat ./testTaj.sh.log
    RVAL=1
fi


echo "Testing SFS"
#./testSFS.sh $WDIR
if [ ! $? -eq 0  ]   ;then
    echo "Problem with SFS exit code: $?"
    cat ./testSFS.sh.log
    RVAL=1
fi

if [[ ! -z "$BAMDIR" ]]; then
echo "Testing basic mpileup"
./testBam.sh $WDIR/angsd $BAMDIR
if [ ! $? -eq 0  ]   ;then
    echo "Problem with basic pileup exit code: $?"
    cat ./testBam.sh.log
    RVAL=1
fi
fi
exit ${RVAL}

if false; then
    echo ./testBam.sh $PRG
    ./testBam.sh $PRG
    
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
