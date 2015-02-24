ANGSD='../angsd0.446/angsd'
echo "==========$0==============="
if [ $# -eq 1 ] 
then
    ANGSD=$1


echo "Using binary: " $ANGSD 
OUTFILE="temp.txt"
echo "output from the program stored in "$OUTFILE

rm -f $OUTFILE
fi

if [ $# -eq 1 ] 
then


rm -f output/fasta*
IND=/home/software/angsd/test/smallBam/smallNA12763.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam

$ANGSD -out output/fasta1 -nThreads 2 -i $IND -doFasta 1 -explode 1 &>>$OUTFILE &
$ANGSD -out output/fasta2 -nThreads 2 -i $IND -doFasta 2 -doCounts 1 -explode 1 &>>$OUTFILE 
$ANGSD -out output/fasta3 -nThreads 2 -i $IND -doFasta 3  -explode 1 &>>$OUTFILE


################
cmp oldResults/fasta1.fa.gz output/fasta1.fa.gz 
cmp oldResults/fasta2.fa.gz output/fasta2.fa.gz 
cmp oldResults/fasta3.fa.gz output/fasta3.fa.gz 

echo "end of list. bang BANG"

#####################################################################
#    make old results
######################################################################
fi
if [ $# -eq 0 ] 
then
IND=/home/software/angsd/test/smallBam/smallNA12763.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam

#ANGSD='../angsd0.563/angsd'

$ANGSD -out oldResults/fasta1 -nThreads 2 -i $IND -doFasta 1 &>>/dev/null &
$ANGSD -out oldResults/fasta2 -nThreads 2 -i $IND -doFasta 2 -doCounts 1 &>>/dev/null &
##bug in -doFasta3 fixed in 0.586
ANGSD='../angsd0.586/angsd'
$ANGSD -out oldResults/fasta3 -nThreads 2 -i $IND -doFasta 3  &>>/dev/null

fi

echo "==========$0==============="
