##modified for 0.576
ANGSD='../angsd0.446/angsd'
echo "==========$0==============="
if [ $# -eq 1 ] 
then
    ANGSD=$1

BEAGLE=/home/albrecht/bin/prog/beagle/beagle.jar

#echo "Using binary: " $ANGSD 
OUTFILE="temp.txt"
#echo "output from the program stored in "$OUTFILE

rm -f $OUTFILE
fi

if [ $# -eq 1 ] 
then

rm -f output/CEU*

$ANGSD -GL 1 -out output/CEU330 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 24 -doMaf 1 -bam CEU330.filelist -r 1:14000000-14010000  -yBin oldResults/CEU330beagle.ybin -doAsso 2 -cov oldResults/CEU330beagle.cov -doPost 1 &>>$OUTFILE

#asso su yeon
$ANGSD -GL 1 -out output/CEU330A1 -nThreads 10 -doMajorMinor 1 -SNP_pval 24 -doMaf 1  -bam CEU330.filelist  -r 1:14000000-14010000 -yBin oldResults/CEU330beagle.ybin -doAsso 1  &>>$OUTFILE

#do imputation
java -Xmx15000m -jar $BEAGLE like=output/CEU330.beagle.gz out=output/CEU330beagle  &>>$OUTFILE

$ANGSD -out output/CEU330beagle -doMaf 4 -beagle output/CEU330beagle.CEU330.beagle.gz.gprobs.gz  -fai hg19.fa.fai &>>$OUTFILE

$ANGSD -out output/CEU330beagleA2  -beagle output/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yQuant oldResults/CEU330beagle.y -doAsso 2 -doMaf 4 -fai hg19.fa.fai &>>$OUTFILE

$ANGSD -out output/CEU330beagleA3  -beagle output/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yBin oldResults/CEU330beagle.ybin -doAsso 2 -doMaf 4  -fai hg19.fa.fai &>>$OUTFILE

$ANGSD -out output/CEU330beagleA4  -beagle output/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yBin oldResults/CEU330beagle.ybin -doAsso 2 -doMaf 4 -cov oldResults/CEU330beagle.cov -fai hg19.fa.fai  &>>$OUTFILE

fi

echo "List of errors:"

function f(){ 
    cut -f2 --complement $1|sed 1d
}

function fz(){ 
    gunzip -c $1 | cut -f2 --complement |sed 1d
}

function f2(){ 
    cut -f2 --complement $1|sed 1d|cut -f1-6
}
function f2z(){ 
    gunzip -c $1 | cut -f2 --complement |sed 1d|cut -f1-6
}

################
cmp output/CEU330beagle.CEU330.beagle.gz.gprobs.gz oldResults/CEU330beagle.CEU330.beagle.gz.gprobs.gz
cmp <(fz output/CEU330beagle.mafs.gz) <(f oldResults/CEU330beagle.mafs)
cmp <(zcat output/CEU330.lrt0.gz) <(zcat oldResults/CEU330.lrt0.gz )
cmp output/CEU330A1.lrt0.gz oldResults/CEU330A1.lrt0.gz
cmp <(f2z output/CEU330beagleA2.lrt0.gz)  <(f2 oldResults/CEU330beagleA2.lrt0 )
cmp <(f2z output/CEU330beagleA3.lrt0.gz) <(f2 oldResults/CEU330beagleA3.lrt0 )
cmp <(zcat output/CEU330beagleA4.lrt0.gz) <(zcat oldResults/CEU330beagleA4.lrt0.gz )

echo "end of list. bang BANG"

#####################################################################
#    make old results
######################################################################

if [ $# -eq 0 ] 
then
BEAGLE=/home/albrecht/bin/prog/beagle/beagle.jar
Rscript -e 'write.table(cbind(rnorm(330)),file="oldResults/CEU330beagle.y",col=F,row=F,qu=F)'

Rscript -e 'write.table(cbind(rbinom(330,1,0.5)),file="oldResults/CEU330beagle.ybin",col=F,row=F,qu=F)'
Rscript -e 'write.table(cbind(rbinom(330,1,0.5),rnorm(330)),file="oldResults/CEU330beagle.cov",col=F,row=F,qu=F)'


for f in `seq 1 10`
do 
cat smallBam.filelist
done > CEU330.filelist


############ new method for logistic regression
#ANGSD='../angsd0.511/angsd'
#$ANGSD -GL 1 -out oldResults/CEU330 -nThreads 10 -doGlf 2 -doMajorMinor 1   -minLRT 24 -doMaf 2 -doSNP 2 -bam CEU330.filelist -r 1:14000000-14010000  -yBin oldResults/CEU330beagle.ybin -doAsso 2 -cov oldResults/CEU330beagle.cov -doPost 1
ANGSD='../angsd0.613/angsd'
 $ANGSD -GL 1 -out oldResults/CEU330 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 24 -doMaf 1 -bam CEU330.filelist -r 1:14000000-14010000  -yBin oldResults/CEU330beagle.ybin -doAsso 2 -cov oldResults/CEU330beagle.cov -doPost 1
$ANGSD -out oldResults/CEU330beagleA4  -beagle output/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yBin oldResults/CEU330beagle.ybin -doAsso 2 -doMaf 4 -cov oldResults/CEU330beagle.cov -fai hg19.fa.fai  


#asso su yeon

#$ANGSD -GL 1 -out oldResults/CEU330A1 -nThreads 10 -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 2 -bam CEU330.filelist  -r 1:14000000-14010000 -yBin oldResults/CEU330beagle.ybin -doAsso 1 
ANGSD='../angsd0.581/angsd'
$ANGSD -GL 1 -out oldResults/CEU330A1 -nThreads 10 -doMajorMinor 1 -SNP_pval 24 -doMaf 1 -bam CEU330.filelist  -r 1:14000000-14010000 -yBin oldResults/CEU330beagle.ybin -doAsso 1 



ANGSD='../angsd0.446/angsd'
#make beagle
#$ANGSD -GL 1 -out oldResults/CEU330 -nThreads 10 -doGlf 2 -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 2 -bam CEU330.filelist -r 1:14000000-14010000  -yBin oldResults/CEU330beagle.ybin -doAsso 2 -cov oldResults/CEU330beagle.cov -doPost 1


#do imputation
java -Xmx15000m -jar $BEAGLE like=oldResults/CEU330.beagle.gz out=oldResults/CEU330beagle

$ANGSD -out oldResults/CEU330beagle -doMaf 16 -beagle oldResults/CEU330beagle.CEU330.beagle.gz.gprobs.gz



#new print format 0123->ACGT
ANGSD='../angsd0.511/angsd'
$ANGSD -out oldResults/CEU330beagleA2  -beagle oldResults/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yQuant oldResults/CEU330beagle.y -doAsso 2 -doMaf 16 -fai hg19.fa.fai

$ANGSD -out oldResults/CEU330beagleA3  -beagle oldResults/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yBin oldResults/CEU330beagle.ybin -doAsso 2 -doMaf 16 -fai hg19.fa.fai


#$ANGSD -out oldResults/CEU330beagleA2  -beagle oldResults/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yQuant oldResults/CEU330beagle.y -doAsso 2 -doMaf 16

#$ANGSD -out oldResults/CEU330beagleA3  -beagle oldResults/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yBin oldResults/CEU330beagle.ybin -doAsso 2 -doMaf 16

#$ANGSD -out oldResults/CEU330beagleA4  -beagle oldResults/CEU330beagle.CEU330.beagle.gz.gprobs.gz -yBin oldResults/CEU330beagle.ybin -doAsso 2 -doMaf 16 -cov oldResults/CEU330beagle.cov



fi

echo "==========$0==============="