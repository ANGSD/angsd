if [[ ! $# -eq 1 ]]
then
    echo "Supply the WORKDIR"
    exit 1
exit
fi

LOG=${0}.log
rm -f ${LOG}

WDIR=$1
ASSODIR=assotest

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 2 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso2new=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 4 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso4new=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 5 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso5new=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yQuant ${ASSODIR}/test.phe -doAsso 6 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso6new=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 2 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso2binnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 4 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso4binnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 5 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso5binnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -beagle ${ASSODIR}/test.beagle -fai hg19.fa.fai  -yBin ${ASSODIR}/test.phe -doAsso 6 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso6binnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -bgen ${ASSODIR}/testB13zlib.bgen -fai hg19.fa.fai  -sampleFile ${ASSODIR}/test.sample -doAsso 2 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso2binbgennew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -bgen ${ASSODIR}/testB13zlib.bgen -fai hg19.fa.fai  -sampleFile ${ASSODIR}/test.sample -doAsso 4 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso4binbgennew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -bgen ${ASSODIR}/testB13zlib.bgen -fai hg19.fa.fai  -sampleFile ${ASSODIR}/test.sample -doAsso 5 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso5binbgennew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -bgen ${ASSODIR}/testB13zlib.bgen -fai hg19.fa.fai  -sampleFile ${ASSODIR}/test.sample -doAsso 6 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso6binbgennew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`


tmp=`$WDIR/angsd -doMaf 4 -vcf-gp ${ASSODIR}/test.vcf -yQuant ${ASSODIR}/test.phe -doAsso 2 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso2vcfnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -vcf-gp ${ASSODIR}/test.vcf -yQuant ${ASSODIR}/test.phe -doAsso 4 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso4vcfnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -vcf-gp ${ASSODIR}/test.vcf -yQuant ${ASSODIR}/test.phe -doAsso 5 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso5vcfnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

tmp=`$WDIR/angsd -doMaf 4 -vcf-gp ${ASSODIR}/test.vcf -yQuant ${ASSODIR}/test.phe -doAsso 6 -out tmp -minCount 0 -minHigh 0 -seed 123 2>>${LOG}`
asso6vcfnew=`zcat tmp.lrt0.gz | tail -n +2 | cut -f7`

## check results LRT and beta for each of the 12 analyses give back RVAL which one that fails

RVAL=0

asso2=1.607704

if [ ! "$asso2" = "$asso2new"  ] ;then
    echo "--------------"
    echo "Problem with score test - beagle file"
    echo "--------------"
    RVAL=2
fi

asso4=1.606574

if [ ! "$asso4" = "$asso4new"  ] ;then
    echo "--------------"
    echo "Problem with latent genotype test - beagle file"
    echo "--------------"
    RVAL=3
fi

asso5=1.607704

if [ ! "$asso5" = "$asso5new"  ] ;then
    echo "--------------"
    echo "Problem with hybrid test - beagle file"
    echo "--------------"
    RVAL=4
fi

asso6=1.615975

if [ ! "$asso6" = "$asso6new"  ] ;then
    echo "--------------"
    echo "Problem with dosage test - beagle file"
    echo "--------------"
    RVAL=5
fi

asso2bin=1.620294

if [ ! "$asso2bin" = "$asso2binnew"  ] ;then
    echo "--------------"
    echo "Problem with binary score test - beagle file"
    echo "--------------"
    RVAL=6
fi

asso4bin=1.540588

if [ ! "$asso4bin" = "$asso4binnew"  ] ;then
    echo "--------------"
    echo "Problem with binary latent genotype test - beagle file"
    echo "--------------"
    RVAL=7
fi

asso5bin=1.620294

if [ ! "$asso5bin" = "$asso5binnew"  ] ;then
    echo "--------------"
    echo "Problem with binary hybrid test - beagle file"
    echo "--------------"
    RVAL=8
fi

asso6bin=1.540588

if [ ! "$asso6bin" = "$asso6binnew"  ] ;then
    echo "--------------"
    echo "Problem with binary dosage test - beagle file"
    echo "--------------"
    RVAL=9
fi

asso2binbgen=1.620294

if [ ! "$asso2binbgen" = "$asso2binbgennew"  ] ;then
    echo "--------------"
    echo "Problem with binary score test - bgen file"
    echo "--------------"
    RVAL=10
fi

asso4binbgen=1.540588

if [ ! "$asso4binbgen" = "$asso4binbgennew"  ] ;then
    echo "--------------"
    echo "Problem with binary latent genotype test - bgen file"
    echo "--------------"
    RVAL=11
fi

asso5binbgen=1.620294

if [ ! "$asso5binbgen" = "$asso5binbgennew"  ] ;then
    echo "--------------"
    echo "Problem with binary hybrid test - bgen file"
    echo "--------------"
    RVAL=12
fi

asso6binbgen=1.540588

if [ ! "$asso6binbgen" = "$asso6binbgennew"  ] ;then
    echo "--------------"
    echo "Problem with binary dosage test - bgen file"
    echo "--------------"
    RVAL=13
fi

asso2=1.607704

if [ ! "$asso2" = "$asso2vcfnew"  ] ;then
    echo "--------------"
    echo "Problem with score test - VCF file"
    echo "--------------"
    RVAL=14
fi

asso4=1.606574

if [ ! "$asso4" = "$asso4vcfnew"  ] ;then
    echo "--------------"
    echo "Problem with latent genotype test - VCF file"
    echo "--------------"
    RVAL=15
fi

asso5=1.607704

if [ ! "$asso5" = "$asso5vcfnew"  ] ;then
    echo "--------------"
    echo "Problem with hybrid test - VCF file"
    echo "--------------"
    RVAL=16
fi

asso6=1.615975

if [ ! "$asso6" = "$asso6vcfnew"  ] ;then
    echo "--------------"
    echo "Problem with dosage test - VCF file"
    echo "--------------"
    RVAL=17
fi


exit $RVAL
