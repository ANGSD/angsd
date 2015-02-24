#!/bin/bash

PRG=""

if [ $# -eq 0 ] 
then
    exit 1;
fi

if [ $# -eq 1 ]
then
    PRG=$1
fi

WDIR=`dirname $PRG`

./testSFS.sh $WDIR || echo "Problem with SFS exit code: $?"

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
    
    echo ./tajTest6.sh
    ./tajTest6.sh $WDIR
    
    echo ./testFasta.sh
    ./testFasta.sh $PRG
    
    echo ./testBaq.sh $PRG
    ./testBaq3.sh $PRG
    
    echo "Netaccess is now deprecated"
    #echo ./testNetAccess.sh $PRG
    #./testNetAccess.sh $PRG
fi
