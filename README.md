## [![Build Status](https://github.com/ANGSD/angsd/actions/workflows/build-tests.yml/badge.svg)](https://github.com/ANGSD/angsd/actions/workflows/build-tests.yml) 



=====
Program for analysing NGS data. 

https://www.popgen.dk/angsd
 
Installation:
=====
1) Using a local folder containing htslib
```
#download htslib
git clone --recurse-submodules https://github.com/samtools/htslib.git;
#download angsd
git clone https://github.com/angsd/angsd.git;

#install htslib
cd htslib
make

#install angsd
cd ../angsd
make HTSSRC=../htslib
```

2) Systemwide installation of htslib

```
git clone https://github.com/angsd/angsd.git;
cd angsd; make HTSSRC=systemwide
```

3) Using htslib submodule

```
git clone https://github.com/angsd/angsd.git;
cd angsd; make
```



Notes
====
* I've switched over to using htslib for parsing single reads (to allow for CRAM reading, while avoid having to write my own CRAM parser). I'm still using my own readpools. Users should therefore also download and install htslib.
* If you are on a mac computer and the compilation process complains about a missnig crybtolib library then do 'make CRYPTOLIB=""'

VCF/BCF allele convention
====
ANGSD's genotype-likelihood and genotype-probability outputs use the
major/minor alleles inferred for each site. When ANGSD writes a VCF or BCF,
the major allele is written as `REF` and the minor allele as `ALT`; these are
not necessarily the reference-genome and alternate alleles. If an assembly
reference is supplied, ANGSD records that allele separately as `INFO/RA`.

The three likelihood/probability values are ordered as major/major,
major/minor and minor/minor. Downstream software must preserve this ordering
and should check allele orientation before combining ANGSD output with a
reference-oriented VCF.

Program has a paper
=====
http://www.biomedcentral.com/1471-2105/15/356/abstract
