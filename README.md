[![Build Status](https://travis-ci.org/ANGSD/angsd.svg?branch=master)](https://travis-ci.org/ANGSD/angsd)
angsd:
=====
Program for analysing NGS data. 

http://www.popgen.dk/angsd
 
Installation:
=====
1) Using a local folder containing htslib

git clone https://github.com/samtools/htslib.git;

git clone https://github.com/angsd/angsd.git;

cd htslib;make;cd ../angsd;make HTSSRC=../htslib

2) Systemwide installation of htslib

git clone https://github.com/angsd/angsd.git;

cd angsd;make

Notes
====
* I've switched over to using htslib for parsing single reads (to allow for CRAM reading, while avoid having to write my own CRAM parser). I'm still using my own readpools. Users should therefore also download and install htslib.

Program has a paper
=====
http://www.biomedcentral.com/1471-2105/15/356/abstract

