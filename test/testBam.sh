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



#SAM=/home/thorfinn/install/samtools/samtools
#${SAM} mpileup -b smallBam.filelist -r 1:14000000-14030000 -q 0 -Q 0 -Bx 2>>temp.txt |cut -f3 --complement >oldResults/mpile
#md5=`md5sum oldResults/mpile|cut -f1 -d" "`
md5=7af69295d04e8b76ebc2073a0d836884
md5new=`$ANGSD -out temp -bam smallBam.filelist -show 1 -r 1:14000000-14030000 -minMapQ 0 -minQ 0 2>>${LOG}|md5sum |cut -f1 -d" "`
RVAL=0
if [ ! "$md5" = "$md5new"  ] ;then
    echo "--------------"
    echo "Problem with first bam pileup test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=2
fi


#${SAM} mpileup -q 30 -b smallBam.filelist -r 1:14000000-14030000 -Q 0 -Bx 2>>temp.txt|cut -f3 --complement   >oldResults/mpile30q
#md5=`md5sum oldResults/mpile30q|cut -f1 -d" "`
md5=a1cd980300100064e0116bc4a470f787

md5new=`$ANGSD -out temp -bam smallBam.filelist -minMapQ 30 -show 1 -r 1:14000000-14030000 2>>${LOG}|md5sum|cut -f1 -d" "` 

if [ ! "$md5" = "$md5new"  ] ;then
    echo "--------------"
    echo "Problem with second bam pileup test"
    echo "--------------"
    cat ${LOG}
    echo "--------------"
    RVAL=3
fi

exit $RVAL