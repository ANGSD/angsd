#!/bin/bash

MD5=md5sum

if [ $# -eq 1 ] 
then
    WDIR=$1
fi


if [[ "$OSTYPE" == "darwin"* ]]; then
    MD5=./md5osx.sh
fi


LOG=${0}.log
rm -rf saf_beagle/output ${LOG}
echo "Cleaning old output dir saf_beagle" &>>${LOG}
mkdir -p saf_beagle/output
ODIR=saf_beagle/output/

REF=saf_beagle/output/ref.fa
BEAG=saf_beagle/input/beagle.txt
TRUTH=saf_beagle/input/true.saf.print

#generate inputs if not already present
if [[ ! -d saf_beagle/input ]]; then
  mkdir -p saf_beagle/input
  
  echo "Generating beagle input" &>>${LOG}
  printf "marker\tallele1\tallele2\tInd0\tInd0\tInd0\tInd1\tInd1\tInd1\tInd2\tInd2\tInd2\n" >${BEAG}
  #should all make symmetric SAF centered around mid-freq bin, with different 'peakiness'
  printf "chr1_1\t0\t1\t0.1\t0.8\t0.1\t0.1\t0.8\t0.1\t0.1\t0.8\t0.1\n" >>${BEAG}
  printf "chr1_2\t0\t1\t0.8\t0.1\t0.1\t0.1\t0.8\t0.1\t0.1\t0.1\t0.8\n" >>${BEAG}
  printf "chr1_3\t0\t1\t1.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t1.0\n" >>${BEAG}
  #should make uniform SAF
  printf "chr1_4\t0\t1\t0.33\t0.33\t0.33\t0.33\t0.33\t0.33\t0.33\t0.33\t0.33\n" >>${BEAG}
  printf "chr1_5\t0\t1\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\n" >>${BEAG}
  #should make mirror-image SAF with mode at lowest/highest freq
  printf "chr1_6\t0\t1\t0.8\t0.1\t0.1\t0.8\t0.1\t0.1\t0.8\t0.1\t0.1\n" >>${BEAG}
  printf "chr1_7\t0\t1\t0.1\t0.1\t0.8\t0.1\t0.1\t0.8\t0.1\t0.1\t0.8\n" >>${BEAG}

  echo "Generating true SAF" &>>${LOG}
  #symmetric SAF centered around mid-freq bin, with different 'peakiness'
  printf "chr1\t1\t-5.345201\t-3.265759\t-1.405563\t0.000000\t-1.405563\t-3.265759\t-5.345201\n" >${TRUTH}
  printf "chr1\t2\t-1.963610\t-0.851204\t-0.851204\t0.000000\t-0.851204\t-0.851204\t-1.963610\n" >>${TRUTH}
  printf "chr1\t3\t-inf\t-inf\t-inf\t0.000000\t-inf\t-inf\t-inf\n" >>${TRUTH}
  #uniform SAF
  printf "chr1\t4\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\n" >>${TRUTH}
  printf "chr1\t5\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\n" >>${TRUTH}
  #mirror-image SAF with mode at lowest/highest freq
  printf "chr1\t6\t0.000000\t-2.079442\t-3.283414\t-4.589666\t-5.362856\t-6.238325\t-6.238325\n" >>${TRUTH}
  printf "chr1\t7\t-6.238325\t-6.238325\t-5.362856\t-4.589666\t-3.283414\t-2.079442\t0.000000\n" >>${TRUTH}

  #may need to `touch ${REF}.fai`
fi

echo "Generating reference fasta" &>>${LOG}
printf ">chr1\nAAAAAAAAAAAAAAA\n" >${REF}
sleep 2 && printf "chr1\t15\t6\t15\t16\n" >${REF}.fai

echo "Generating saf likelihoods" &>>${LOG}
$WDIR/angsd -beagle ${BEAG} -ref ${REF} -anc ${REF} -fai ${REF}.fai -out ${ODIR}/out -dosaf 4 &>>${LOG}
$WDIR/misc/realSFS print ${ODIR}/out.saf.idx 1>${ODIR}/out.saf.print 2>>${LOG}

cmp ${TRUTH} ${ODIR}/out.saf.print >>${LOG} >&1 || exit 1
