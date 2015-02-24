#!/bin/bash


echo "==========$0==============="
echo "You should supply the rootdir of the sourcepackage, this analysis should take around 2 minutes"

WDIR=../../angsd0.556/

if [ $# -eq 1 ] 
then
    WDIR=$1
fi

echo "Using folder: " $WDIR

OS="temp.txt"
echo "Cleaning old output dir 'taj/'"
rm -rf taj
mkdir -p taj
echo "Generating simulated haplotypes"
msms -ms 40 100 -t 900 -r 400 -SAA 1000 -SaA 500 -N 10000 -SF 0 -Sp .5 -oTPi -oAFS >taj/msout
grep -i afs taj/msout |grep -i sum|cut -f2 -d":"|tr " " "\n" >taj/msout.sfs
echo "Generating genotype likelihood based on the haplotypes"
${WDIR}/misc/msToGlf -in taj/msout -out taj/glout -err 0.005 -depth 8 -singleOut 1 2>${OS}
echo "1) Will do unfolded analysis"
${WDIR}/angsd -isSim 1 -glf taj/glout.glf.gz -out taj/norm -doSaf 1 -nInd 20 -fai hg19.fa.fai  2>${OS}
${WDIR}/misc/emOptim2 taj/norm.saf 40 -P 24 -nSites 1000000 >taj/norm.saf.em.ml 2>${OS}
${WDIR}/angsd -isSim 1 -glf taj/glout.glf.gz -out taj/norm -nInd 20 -doThetas 1 -doSaf 1 -pest taj/norm.saf.em.ml -fai hg19.fa.fai  2>${OS}
${WDIR}/misc/thetaStat make_bed taj/norm.thetas.gz  2>${OS}
${WDIR}/misc/thetaStat do_stat taj/norm.thetas.gz -nChr 40  2>${OS}

echo "2) Will do folded analysis"

${WDIR}/angsd -isSim 1 -glf taj/glout.glf.gz -out taj/fold -doSaf 1 -nInd 20 -fold 1 -fai hg19.fa.fai  2>${OS}
${WDIR}/misc/emOptim2 taj/fold.saf 20 -P 24 -nSites 1000000 >taj/fold.saf.em.ml 2>${OS}
${WDIR}/angsd -isSim 1 -glf taj/glout.glf.gz -out taj/fold -nInd 20 -doThetas 1 -doSaf 1 -pest taj/fold.saf.em.ml -fold 1  -fai hg19.fa.fai 2>${OS}
${WDIR}/misc/thetaStat make_bed taj/fold.thetas.gz  2>${OS}
${WDIR}/misc/thetaStat do_stat taj/fold.thetas.gz -nChr 40  2>${OS}

echo -e "\tThe theta estimates from msms simulated haplotypes for 8x 0.5% errorrate"
echo -e "\tWatterson\tpairwise\ttajimasD"
grep ThetaW taj/msout|grep Sum
echo -e "-------\nThe estimated thetas (sum of 100 reps) and the tajima value (unfolded)"
cat taj/norm.thetas.gz.pestPG
echo -e "-------\nThe estimated thetas (sum of 100 reps) and the tajima value (unfolded)"
cat taj/fold.thetas.gz.pestPG
echo "You should eyeball the above and see if they are comparable column (1-5,9) (not all defined for unfold)"


echo "==========$0==============="


