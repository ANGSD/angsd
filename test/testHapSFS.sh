#!/bin/bash
MD5=md5sum

if [ $# -eq 1 ] 
then
    WDIR=$1
fi

TDIR=haploid_SFS

LOG=${0}.log
echo "Cleaning old output dir ${TDIR}/" &>${LOG}
rm -rf ${TDIR}/output ${LOG}
mkdir -p ${TDIR}/output

MSMS=${TDIR}/input/msout.gz
gunzip -c ${MSMS} >${TDIR}/input/msout
MSMS=${TDIR}/input/msout
ODIR=${TDIR}/output/
##SIMfiles from tajima test

echo "Generating genotype likelihood based on the haplotypes" >>${LOG} 2>&1
${WDIR}/misc/msToGlf -in ${MSMS} -out ${ODIR}/glout -err 0.005 -depth 8 -singleOut 1 -simpleRand 0 -simhap 1 >>${LOG} 2>&1
${WDIR}/angsd -isSim 1 -ishap 1 -glf ${ODIR}/glout.glf.gz -out ${ODIR}/norm -doSaf 1 -nInd 40 -fai hg19.fa.fai  2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/norm.saf.idx -P 4  -m 0 -seed -1 >${ODIR}/norm.saf.em.ml 2>>${LOG}
${WDIR}/misc/realSFS saf2theta ${ODIR}/norm.saf.idx -sfs ${ODIR}/norm.saf.em.ml -outname ${ODIR}/norm.thetaFromSaf 2>>${LOG}

${WDIR}/misc/thetaStat do_stat ${ODIR}/norm.thetas.idx 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/norm.thetas.idx -outnames ${ODIR}/norm.thetas.idx.win -step 1000 -win 5000 2>>${LOG}

${WDIR}/misc/thetaStat do_stat ${ODIR}/norm.thetaFromSaf.thetas.idx 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/norm.thetaFromSaf.thetas.idx -outnames ${ODIR}/norm.thetaFromSaf.thetas.idx.win -step 1000 -win 5000 2>>${LOG}

echo "2) Will do folded analysis" >>${LOG} 2>&1
${WDIR}/misc/realSFS ${ODIR}/norm.saf.idx -P 4  -m 0 -seed -1 -fold 1 >${ODIR}/fold2.saf.em.ml 2>>${LOG}
${WDIR}/misc/realSFS saf2theta ${ODIR}/norm.saf.idx -sfs ${ODIR}/fold2.saf.em.ml -outname ${ODIR}/fold.thetaFromSaf -fold 1 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/fold.thetaFromSaf.thetas.idx -nChr 40 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/fold.thetaFromSaf.thetas.idx -nChr 40 -outnames ${ODIR}/fold.thetaFromSaf.thetas.idx.win -step 1000 -win 5000 2>>${LOG}


echo -e "\tThe theta estimates from msms simulated haplotypes for 8x 0.5% errorrate"  >>${LOG} 2>&1
echo -e "\tWatterson\tpairwise\ttajimasD" >>${LOG} 2>&1
grep ThetaW ${MSMS}|grep Sum >>${LOG} 2>&1
echo -e "-------\nThe estimated thetas (sum of 100 reps) and the tajima value (unfolded)" >>${LOG} 2>&1
cat ${ODIR}/norm.thetaFromSaf.thetas.idx.pestPG >>${LOG} 2>&1
echo -e "-------\nThe estimated thetas (sum of 100 reps) and the tajima value (folded)" >>${LOG} 2>&1
cat ${ODIR}/fold.thetaFromSaf.thetas.idx.pestPG >>${LOG} 2>&1
echo "You should eyeball the above and see if they are comparable column (1-5,9) (not all defined for unfold)" >>${LOG} 2>&1

##when generated the results:
##md5sum tajima/output/*pestPG tajima/output/*.ml tajima/output/*.saf tajima/output/*.gz  >tajima/md5/pestPG.md5sum
#${MD5}  -c tajima/md5/pestPG.md5sum >>${LOG} 2>&1 || exit 1
#${MD5}  -c ${TDIR}/md5/pestPG3.md5sum >>${LOG} 2>&1 || exit 2
 
