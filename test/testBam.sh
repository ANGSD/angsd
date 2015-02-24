echo "==========$0=(seconds)=============="
if [[ $# == 0 ]]
then
    echo "Supply an angsd binary"
exit
fi

rm -f oldResults/mpile output/mpile oldResults/mpile30q output/mpile30q

ANGSD=$1
SAM=/opt/samtools-0.1.19/samtools

${SAM} mpileup -b smallBam.filelist -r 1:14000000-14030000 -q 0 -Q 0 2>>temp.txt |cut -f3 --complement >oldResults/mpile
$ANGSD -out temp -bam smallBam.filelist -show 1 -r 1:14000000-14030000  >output/mpile 2>>temp.txt


${SAM} mpileup -q 30 -b smallBam.filelist -r 1:14000000-14030000 -Q 0 2>>temp.txt|cut -f3 --complement   >oldResults/mpile30q
$ANGSD -out temp -bam smallBam.filelist -minMapQ 30 -show 1 -r 1:14000000-14030000 2>>temp.txt  >output/mpile30q

echo "List of errors:"
#dont print N for the reference if not supplied
cmp oldResults/mpile output/mpile
cmp oldResults/mpile30q output/mpile30q
#md5sum oldResults/mpile output/mpile oldResults/mpile30q output/mpile30q
echo "end of list. bang BANG"

echo "==========$0==============="
