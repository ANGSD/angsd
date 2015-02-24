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


rm -f output/abba*

$ANGSD -out output/abba -nThreads 1 -bam smallerBam.filelist -doAbbababa 1 -doCounts 1  -anc /space/genomes/refgenomes/ancestral/hg19/fasta/hg19ancNoChr.fa -r 1: &>>$OUTFILE


################
cmp output/abba.abbababa oldResults/abba.abbababa

echo "end of list. bang BANG"

#####################################################################
#    make old results
######################################################################
fi
if [ $# -eq 0 ] 
then

ANGSD='../angsd0.562/angsd'
$ANGSD -out oldResults/abba -nThreads 1 -bam smallerBam.filelist -doAbbababa 1 -doCounts 1  -anc /space/genomes/refgenomes/ancestral/hg19/fasta/hg19ancNoChr.fa  -r 1: &>/dev/null 

fi

echo "==========$0==============="