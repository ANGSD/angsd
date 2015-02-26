#!/bin/bash
echo "==========$0=(10minutes)=============="
PRG=""
if [[ $# == 0 ]]
then
    echo "Supply an angsd binary"
exit
fi

PRG=$1
PRG2="/opt/samtools-0.1.19/samtools"
FA="/space/genomes/refgenomes/hg19/tsk/hg19.fa"


if [ ! -f ${PRG} ];
then
    echo "Angsd binary not found: \"${PRG}\" will exit"
    exit;
fi



echo "Using binary: " ${PRG} 
echo  "TEST 1/3 : Checking no reference, no baq "
rm -f output/sam1.mpil output/an1.mpil
${PRG2} mpileup -b  smallBam.filelist -q 0 -Q 0 2>>temp.txt|cut -f3 --complement >output/sam1.mpil
${PRG} -out here -bam smallBam.filelist -show 1 >output/an1.mpil 2>>temp.txt

cmp output/sam1.mpil output/an1.mpil
md5sum output/sam1.mpil output/an1.mpil
echo "Done checking first run"
#exit

echo "TEST 2/3 Checking with reference, no baq"
rm -f output/sam1.mpil output/an1.mpil
${PRG2} mpileup -b  smallBam.filelist -q 0 -Q 0 -f ${FA} -B 2>>temp.txt >output/sam1.mpil
${PRG} -out here -bam smallBam.filelist -show 1 -ref ${FA} -baq 0 2>>temp.txt >output/an1.mpil


cmp output/sam1.mpil output/an1.mpil
md5sum output/sam1.mpil output/an1.mpil
echo "Done checking second run"






echo "TEST 3/3 with reference, normal baq;"
##samtools has problems in the begenning of the reads. There we do active lookup with samtools
rm -f output/sam1.mpil output/an1.mpil
${PRG} -out here -bam smallBam.filelist -show 1 -ref ${FA} -baq 1 >output/an1.mpil 2>>temp.txt
touch output/sam1.mpil

for i in `cut -f1 output/an1.mpil |uniq`
do 
    #echo "Building output (${i}/20)"
    ${PRG2} mpileup -b smallBam.filelist  -f ${FA}   -q 0 -Q 0 -r ${i}: >>output/sam1.mpil 2>temp.txt
done



cmp output/sam1.mpil output/an1.mpil
md5sum output/sam1.mpil output/an1.mpil
echo "Done checking last run"
echo "==========$0==============="