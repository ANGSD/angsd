##modified for 0.576+0.584 (number of effective sites)
echo "==========$0==============="

mangeANGSD='../angsd0.445/angsd'
if [ $# -eq 1 ] 
then
    ANGSD=$1
#    echo "Using binary: " $ANGSD 
    OUTFILE="temp.txt"
#    echo "output from the program stored in "$OUTFILE
fi

if [ $# -eq 1 ] 
then
    rm -f output/smallBamMics*
    
    rm -f oldResults/*.bin oldResults/*.idx

    ANC=/pontus/genomes/refgenomes/ancestral/hg19/fasta/hg19ancNoChr.fa
    REF=/pontus/genomes/refgenomes/perfect/fasta/hg19perfectQ20q30noChr.fa
    ######################
    #ancestral errors
    $ANGSD -doAncError 1 -anc $ANC -ref $REF -out output/smallBamMicsAnc -bam smallBam.filelist -r 1: -minQ 0 &>>$OUTFILE
    #Rscript ../angsd0.41/R/estError.R file=output/smallBamMicsAnc.ancError out=output/smallBamMicsAnc

    $ANGSD -doAncError 1 -anc ${ANC} -ref ${REF} -out output/smallBamMicsAncq30Q20 -bam smallBam.filelist -r 1: -minQ 20 -minMapQ 30 &>>$OUTFILE
    #Rscript ../angsd0.41/R/estError.R file=output/smallBamMicsAnc.ancError out=output/smallBamMicsAncq30Q20.ancError
    ######################
    #beagle 
    $ANGSD -GL 1 -out output/smallBamMicsBeagle -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doGlf 2 -r 1:14000000-14010000 -SNP_pval 24 &>>$OUTFILE
    java -Xmx1000m -jar /home/albrecht/bin/prog/beagle/beagle.jar like=output/smallBamMicsBeagle.beagle.gz out=output/smallBamMicsBeagleOut &>>$OUTFILE
    $ANGSD -beagle output/smallBamMicsBeagleOut.smallBamMicsBeagle.beagle.gz.gprobs.gz -out output/smallBamMicsBeagleOut.smallBamMicsBeagle -doMaf 4 -fai  hg19.fa.fai &>>$OUTFILE


#########################
# filters
#LRT 
    $ANGSD -GL 1 -out output/smallBamMics1 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -r 1: -SNP_pval 14 -doCounts 1 -dumpCounts 3 -minQ 0 -minMapQ 0 &>>$OUTFILE

#quality scores
    $ANGSD -GL 1 -out output/smallBamMics2 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -r 1:14000000-14010000 -doCounts 1 -dumpCounts 4 -minQ 20 -minMapQ 30 &>>$OUTFILE

# -filter
    rm -f oldResults/smallBamMics1.filt2.*
    $ANGSD sites index oldResults/smallBamMics1.filt2
    $ANGSD -GL 1 -r 1: -out output/smallBamMics3 -nThreads 10 -doMajorMinor 1 -doMaf 1 -bam smallBam.filelist -doCounts 1 -dumpCounts 3 -sites oldResults/smallBamMics1.filt2 &>>$OUTFILE


#######################
#test change of chromosomes work
    $ANGSD -bam smallBam.filelist -nInd 4 -doCounts 1  -dumpCounts 1 -out output/tsk.new &>>$OUTFILE 




################
echo "List of errors:"

#ancestral
#echo "ancestral"
cmp output/smallBamMicsAnc.ancError oldResults/smallBamMicsAns.ancError
cmp output/smallBamMicsAncq30Q20.ancError oldResults/smallBamMicsAnsq30Q20.ancError                                                  
 
#md5sum output/smallBamMicsAnc.ancError oldResults/smallBamMicsAns.ancError
#md5sum output/smallBamMicsAncq30Q20.ancError oldResults/smallBamMicsAnsq30Q20.ancError

#beagle
#echo "beagle"
cmp oldResults/smallBamMicsBeagleOut.smallBamMicsBeagle.mafs.gz output/smallBamMicsBeagleOut.smallBamMicsBeagle.mafs.gz
#md5sum oldResults/smallBamMicsBeagleOut.smallBamMicsBeagle.mafs.gz output/smallBamMicsBeagleOut.smallBamMicsBeagle.mafs.gz

#filter
#echo "counts and mafs"
#echo "####counts1"
cmp oldResults/smallBamMics1.counts.gz output/smallBamMics1.counts.gz
cmp oldResults/smallBamMics1.pos.gz output/smallBamMics1.pos.gz
cmp <(gunzip -c oldResults/smallBamMics1.mafs.gz|cut -f6 --complement) <(gunzip -c  output/smallBamMics1.mafs.gz|cut -f6 --complement)

#echo "####counts2"
cmp oldResults/smallBamMics2.counts.gz output/smallBamMics2.counts.gz
cmp oldResults/smallBamMics2.pos.gz output/smallBamMics2.pos.gz
cmp oldResults/smallBamMics2.mafs.gz output/smallBamMics2.mafs.gz

#echo "####counts3"
cmp oldResults/smallBamMics3.counts.gz output/smallBamMics3.counts.gz
cmp oldResults/smallBamMics3.pos.gz output/smallBamMics3.pos.gz
cmp oldResults/smallBamMics3.mafs.gz output/smallBamMics3.mafs.gz

#echo -e "####checking chr changes\n"
cmp oldResults/tsk.old.pos.gz output/tsk.new.pos.gz


#echo "==========$0==============="


echo "end of list. bang BANg"

fi

#####################################################################
#    make old results
######################################################################

if [ $# -eq 0 ] 
then
    rm -f oldResults/smallBamM*

    ANGSD='../angsd0.584/angsd'

    $ANGSD -bam smallBam.filelist -nInd 4 -doCounts 1  -dumpCounts 1 -out oldResults/tsk.old 

    ANGSD='../angsd0.580/angsd'
#ancestral errors
  $ANGSD -doAncError 1 -anc /space/genomes/refgenomes/ancestral/hg19/fasta/hg19ancNoChr.fa -ref /space/genomes/refgenomes/perfect/fasta/hg19perfectQ20q30noChr.fa -out oldResults/smallBamMicsAns -bam smallBam.filelist -r 1: -minQ 0 -sample 0

    $ANGSD -doAncError 1 -anc /space/genomes/refgenomes/ancestral/hg19/fasta/hg19ancNoChr.fa -ref /space/genomes/refgenomes/perfect/fasta/hg19perfectQ20q30noChr.fa -out oldResults/smallBamMicsAnsq30Q20 -bam smallBam.filelist -r 1: -minQ 20 -minMapQ 30 -sample 0

#    $ANGSD -doAnsError 2 -anc /space/genomes/refgenomes/ancestral/hg19/fasta/hg19ancNoChr.fa -ref /space/genomes/refgenomes/perfect/fasta/hg19perfectQ20q30noChr.fa -out oldResults/smallBamMicsAns -bam smallBam.filelist -r 1:
#Rscript ../angsd0.41/R/estError.R file=oldResults/smallBamMicsAns.ancError out=oldResults/smallBamMicsAns

 #   $ANGSD -doAnsError 2 -anc /space/genomes/refgenomes/ancestral/hg19/fasta/hg19ancNoChr.fa -ref /space/genomes/refgenomes/perfect/fasta/hg19perfectQ20q30noChr.fa -out oldResults/smallBamMicsAnsq30Q20 -bam smallBam.filelist -r 1: -minQ 20 -minMapQ 30
#Rscript ../angsd0.41/R/estError.R file=oldResults/smallBamMicsAns.ancError out=oldResults/smallBamMicsAnsq30Q20.ancError

    ANGSD='../angsd0.544/angsd'
#beagle 
    $ANGSD -GL 1 -out oldResults/smallBamMicsBeagle -nThreads 10 -doMajorMinor 1 -doMaf 2 -bam smallBam.filelist -doGlf 2 -r 1:14000000-14010000 -doSNP 1 -minLRT 24
    java -Xmx1000m -jar /home/albrecht/bin/prog/beagle/beagle.jar like=oldResults/smallBamMicsBeagle.beagle.gz out=oldResults/smallBamMicsBeagleOut
    $ANGSD -beagle oldResults/smallBamMicsBeagleOut.smallBamMicsBeagle.beagle.gz.gprobs.gz -out oldResults/smallBamMicsBeagleOut.smallBamMicsBeagle -doMaf 16 -fai  hg19.fa.fai

#filter init
    $ANGSD -GL 1 -out oldResults/smallBamMics1 -nThreads 10 -doMajorMinor 1 -doMaf 2 -bam smallBam.filelist -r 1: -doSNP 1 -minLRT 14 -doCounts 1 -dumpCounts 3 -minQ 0 -minMapQ 0 
    gunzip -c oldResults/smallBamMics1.mafs.gz | sed 1d|cut -f1-2 > oldResults/smallBamMics1.filt

    $ANGSD -GL 1 -out oldResults/smallBamMics2 -nThreads 10 -doMajorMinor 1 -doMaf 2 -bam smallBam.filelist -r 1:14000000-14010000 -minQ 20 -doCounts 1 -dumpCounts 4 -minMapQ 30

# -filter use
    $ANGSD -GL 1 -r 1: -out oldResults/smallBamMics3 -nThreads 10 -doMajorMinor 1 -doMaf 2 -bam smallBam.filelist -doCounts 1 -dumpCounts 3 -filter oldResults/smallBamMics1.filt

 

fi
