#!/bin/bash
# -r/-rf test on one bam from smallBam (indexed, reads at 14000000-14100000 on each chromosome).
# regions/regions.rf has 5 regions on 4 chromosomes, two on the same chromosome, a
# whole chromosome, a comment and an empty line.
# -rf must give the same sites as running each region with -r, in the same order,
# with no duplicated sites, only sites inside the regions, and the same with threads.
# regions/regions.bed has the same regions as BED, regions/mixed.rf as both, and must give
# the same. Bad region files must give exit code 1.
if [[ ! $# -eq 2 ]]
then
	echo "Must supply an angsd binary and the folder with smallBam; will exit"
	exit 1
fi
ANGSD=${1}
BDIR=${2}

LOG=${0}.log
rm -f ${LOG}

echo "==========${0}==========" >> ${LOG}
echo "Running test: ${0} with ${ANGSD}" >> ${LOG}

TDIR=regions
RESDIR=${TDIR}/test_results
rm -rf ${RESDIR}
mkdir -p ${RESDIR}

ls ${BDIR}/smallNA06985.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam > ${RESDIR}/bam.filelist
RF=${TDIR}/regions.rf
ARGS="-bam ${RESDIR}/bam.filelist -doCounts 1 -dumpCounts 1"

RVAL=0
fail(){
	echo "--------------"
	echo "Problem with regions test: $1"
	echo "--------------"
	cat ${LOG}
	echo "--------------"
	RVAL=1
}

#sites of a .pos.gz without the header
body(){
	gunzip -c $1 | tail -n +2
}

#each region by itself with -r
grep -v '^#' ${RF} | grep -v '^$' > ${RESDIR}/regions.txt
rm -f ${RESDIR}/each.pos
i=0
while read reg; do
	i=$((i+1))
	${ANGSD} ${ARGS} -r ${reg} -out ${RESDIR}/r${i} >> ${LOG} 2>&1 || fail "-r ${reg} exit code $?"
	n=$(body ${RESDIR}/r${i}.pos.gz | wc -l)
	if [[ ${n} -eq 0 ]]; then
		fail "no sites in region ${reg}"
	fi
	body ${RESDIR}/r${i}.pos.gz >> ${RESDIR}/each.pos
done < ${RESDIR}/regions.txt

#all regions with -rf, with 1 and 3 threads
for t in 1 3; do
	${ANGSD} ${ARGS} -rf ${RF} -nThreads ${t} -out ${RESDIR}/rf${t} >> ${LOG} 2>&1 || fail "-rf -nThreads ${t} exit code $?"
	body ${RESDIR}/rf${t}.pos.gz > ${RESDIR}/rf${t}.pos
	if ! cmp ${RESDIR}/rf${t}.pos ${RESDIR}/each.pos >> ${LOG} 2>&1; then
		fail "-rf -nThreads ${t} differs from running each region with -r"
	fi
done

#the same regions as BED, gzipped BED, and mixed with chr:start-stop
gzip -c ${TDIR}/regions.bed > ${RESDIR}/regions.bed.gz
for f in ${TDIR}/regions.bed ${RESDIR}/regions.bed.gz ${TDIR}/mixed.rf; do
	n=$(basename ${f})
	${ANGSD} ${ARGS} -rf ${f} -out ${RESDIR}/${n} >> ${LOG} 2>&1 || fail "-rf ${n} exit code $?"
	body ${RESDIR}/${n}.pos.gz > ${RESDIR}/${n}.pos
	if ! cmp ${RESDIR}/${n}.pos ${RESDIR}/each.pos >> ${LOG} 2>&1; then
		fail "-rf ${n} differs from running each region with -r"
	fi
done

#bad region files must fail with exit code 1
printf '1\t14000000\n' > ${RESDIR}/bad_2cols.rf
printf '1\t14000000\t14020000\nchrNotThere\t0\t100\n' > ${RESDIR}/bad_chr_bed.rf
printf 'chrNotThere:1-100\n' > ${RESDIR}/bad_chr.rf
printf '1\t14020000\t14000000\n' > ${RESDIR}/bad_order.rf
printf '1\t14000000\t300000000\n' > ${RESDIR}/bad_end.rf
printf '1\t14000000\tabc\n' > ${RESDIR}/bad_number.rf
printf '# only a comment\n\n' > ${RESDIR}/bad_empty.rf
for f in ${RESDIR}/bad_*.rf; do
	${ANGSD} ${ARGS} -rf ${f} -out ${RESDIR}/bad >> ${LOG} 2>&1
	e=$?
	if [[ ${e} -ne 1 ]]; then
		fail "$(basename ${f}) gave exit code ${e}, not 1"
	fi
done
${ANGSD} ${ARGS} -r chrNotThere -out ${RESDIR}/bad >> ${LOG} 2>&1
e=$?
if [[ ${e} -ne 1 ]]; then
	fail "-r chrNotThere gave exit code ${e}, not 1"
fi

#overlapping regions (also with -r) are used as they are, but must give a warning with the number
#of overlapping regions. 1:14025001-14026000 only touches the region before it, which is not an overlap
printf '1:14000001-14010000\n3:14000001-14005000\n1:14005001-14020000\n1\t14015000\t14025000\n1\t14025000\t14026000\n' > ${RESDIR}/overlap.rf
printf '1:14000001-14025000\n3:14000001-14006000\n1:14025001-14026000\n' > ${RESDIR}/touch.rf
${ANGSD} ${ARGS} -rf ${RESDIR}/overlap.rf -r 3:14004001-14006000 -out ${RESDIR}/overlap > ${RESDIR}/overlap.log 2>&1 || fail "-rf overlap.rf exit code $?"
${ANGSD} ${ARGS} -rf ${RESDIR}/touch.rf -out ${RESDIR}/touch > ${RESDIR}/touch.log 2>&1 || fail "-rf touch.rf exit code $?"
cat ${RESDIR}/overlap.log ${RESDIR}/touch.log >> ${LOG}
if ! grep -q "Warning: 3 regions overlap other regions" ${RESDIR}/overlap.log; then
	fail "no warning (or a wrong number) for 3 overlapping regions"
fi
if grep -q "regions overlap" ${RESDIR}/touch.log; then
	fail "warning for regions that don't overlap"
fi

#no site twice
ndup=$(cut -f1,2 ${RESDIR}/rf1.pos | sort | uniq -d | wc -l)
if [[ ${ndup} -ne 0 ]]; then
	fail "${ndup} sites are in the -rf output more than once"
fi

#every site is inside one of the regions (chr or chr:start-stop, 1-based)
nout=$(awk 'NR==FNR{split($1,a,/[:-]/); n++; c[n]=a[1]; s[n]=(a[2]==""?1:a[2]); e[n]=(a[3]==""?1e12:a[3]); next}
	{ok=0; for(i=1;i<=n;i++) if($1==c[i]&&$2>=s[i]&&$2<=e[i]) ok=1; if(!ok) bad++}
	END{print bad+0}' ${RESDIR}/regions.txt ${RESDIR}/rf1.pos)
if [[ ${nout} -ne 0 ]]; then
	fail "${nout} sites are outside the regions"
fi

echo "==========$0===============" >> ${LOG}

exit $RVAL
