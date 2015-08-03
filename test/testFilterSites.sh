if [[ ! $# -eq 2 ]]
then
    echo "Supply an angsd binary and BAMdir"
    exit 1
exit
fi

LOG=${0}.log
rm -f ${LOG}

ANGSD=$1
BDIR=$2

echo $ANGSD $BDIR

ls $BDIR/*.bam >smallBam.filelist

$ANGSD sites index sites/s1 2>>${LOG}
$ANGSD -out sites/run1 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 1 -sites sites/s1 2>>${LOG}

$ANGSD sites index sites/s2 -compl 1 2>>${LOG}
$ANGSD -out sites/run2 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 1 -sites sites/s2 2>>${LOG}

$ANGSD sites index sites/s3  2>>${LOG}
$ANGSD -out sites/run3 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 1 -sites sites/s3 2>>${LOG}

$ANGSD sites index sites/s4  2>>${LOG}
$ANGSD -out sites/run4 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 3 -sites sites/s4 2>>${LOG}

$ANGSD sites index sites/s5  2>>${LOG}
$ANGSD -out sites/run5 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 3 -sites sites/s5 2>>${LOG}

$ANGSD sites index sites/s5  2>>${LOG}
$ANGSD -out sites/run6 -bam smallBam.filelist -nind 3 -minMapQ 30 -minQ 20  -gl 1 -domaf 1 -doglf 2 -domajorminor 3 -sites sites/s5 -dovcf 1 -dopost 1 2>>${LOG}


cd sites
md5sum  -c md5orig.orig &>>../${LOG}||exit 8

