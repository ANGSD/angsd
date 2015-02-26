#change from how minor is estimated from counts in 0.584 so this is now updated

ANGSD='../angsd0.584/angsd'
echo "==========$0==============="
if [ $# -eq 1 ] 
then
    ANGSD=$1
#    echo "Using binary: " $ANGSD 
    OUTFILE="temp.txt"
#    echo "output from the program stored in "$OUTFILE
fi

if [ $# -eq 1 ] 
then

rm -f output/smallBamGL*
rm -f -r angsd_tmpdir 

#samtools
$ANGSD -GL 1 -out output/smallBamGL1 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 1 -r 1:14000000-14010000 -nlines 5000  &>>$OUTFILE
#GATK
$ANGSD -GL 2 -out output/smallBamGL2 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 1 -r 1:14000000-14010000  -nlines 5000 &>>$OUTFILE
#do Calibration SOAP
$ANGSD -GL 3 -out output/smallBamGL3 -nThreads 10 -bam smallBam.filelist -r 1:14000000-14010000 -ref /pontus/genomes/refgenomes/hg19/tsk/hg19.db135.removed.noChr.fa  -nlines 5000 -minQ 0 &>>$OUTFILE
#SOAP
$ANGSD -GL 3 -out output/smallBamGL3 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 1 -r 1:14000000-14010000   -nlines 5000  &>>$OUTFILE
#Kim 
$ANGSD -GL 4 -out output/smallBamGL4 -nThreads 10 -doMajorMinor 2 -doMaf 1 -bam smallBam.filelist -doGlf 1 -r 1:14000000-14010000 -doCounts 1 -minQ 20 -minMapQ 13   -nlines 5000  &>>$OUTFILE

# other glf output
$ANGSD -GL 1 -out output/smallBamGL1_2 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 2 -r 1:14000000-14010000  -nlines 5000  &>>$OUTFILE
$ANGSD -GL 1 -out output/smallBamGL1_3 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 3 -r 1:14000000-14010000  &>>$OUTFILE
$ANGSD -GL 1 -out output/smallBamGL1_4 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 4 -r 1:14000000-14010000   -nlines 5000 &>>$OUTFILE




################
echo "List of errors:"
cmp output/smallBamGL1.mafs.gz oldResults/smallBamGL1.mafs.gz
cmp output/smallBamGL2.mafs.gz oldResults/smallBamGL2.mafs.gz
cmp output/smallBamGL3.mafs.gz oldResults/smallBamGL3.mafs.gz
cmp output/smallBamGL4.mafs.gz oldResults/smallBamGL4.mafs.gz

cmp output/smallBamGL1.glf.gz oldResults/smallBamGL1.glf.gz
cmp output/smallBamGL2.glf.gz oldResults/smallBamGL2.glf.gz
cmp output/smallBamGL3.glf.gz oldResults/smallBamGL3.glf.gz
cmp output/smallBamGL4.glf.gz oldResults/smallBamGL4.glf.gz

cmp output/smallBamGL1_2.beagle.gz oldResults/smallBamGL1_2.beagle.gz
cmp output/smallBamGL1_3.glf.gz oldResults/smallBamGL1_3.glf.gz
cmp output/smallBamGL1_4.glf.gz oldResults/smallBamGL1_4.glf.gz
echo "end of list. bang BANG"
fi
#####################################################################
#    make old results
######################################################################

if [ $# -eq 0 ] 
then
rm oldResults/smallBamGL1*

ANGSD='../angsd0.584/angsd'
$ANGSD -GL 1 -out oldResults/smallBamGL1 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 1 -r 1:14000000-14010000 
#GATK
$ANGSD -GL 2 -out oldResults/smallBamGL2 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 1 -r 1:14000000-14010000 
#do Calibration SOAP
rm -r angsd_tmpdir 
$ANGSD -GL 3 -out oldResults/smallBamGL1 -nThreads 10 -bam smallBam.filelist -r 1:14000000-14010000 -ref /pontus/genomes/refgenomes/hg19/tsk/hg19.db135.removed.noChr.fa -nlines 5000 -minQ 0
$ANGSD -GL 3 -out oldResults/smallBamGL3 -nThreads 1 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 1 -r 1:14000000-14010000  

$ANGSD -GL 4 -out oldResults/smallBamGL4 -nThreads 10 -bam smallBam.filelist -r 1:14000000-14010000 -doCounts 1 -minQ 20 -minMapQ 13 -doMajorMinor 2 -doMaf 1 -doGlf 1
# other glf output
$ANGSD -GL 1 -out oldResults/smallBamGL1_2 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 2 -r 1:14000000-14010000
$ANGSD -GL 1 -out oldResults/smallBamGL1_3 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 3 -r 1:14000000-14010000
$ANGSD -GL 1 -out oldResults/smallBamGL1_4 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 4 -r 1:14000000-14010000

fi 

echo "==========$0==============="