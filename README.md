angsd
=====
Program for analysing NGS data. 

=====
I've switched over to using htslib for parsing single reads (to allow
for CRAM reading, while avoid having to write my own CRAM parser)
I'm still using my own readpools. Users should therefore also download
and install htslib.

=====

git clone https://github.com/samtools/htslib.git;git clone https://github.com/angsd/angsd.git
cd angsd;make HTSDIR=../htslib;

=====
Manual for ANGSD 
http://www.popgen.dk/angsd/


=====
Program has a paper
http://www.biomedcentral.com/1471-2105/15/356/abstract

