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

Reference and ancestral-allele requirements
=====
ANGSD modes do not all use the same allele convention or require the same external sequence.

* `-anc` supplies an ancestral sequence and is required for polarized spectrum calculations. Without it, use an explicitly folded spectrum when that is scientifically appropriate.
* `-ref` supplies a reference sequence for read alignment and reference-coordinate operations. It does not by itself identify the ancestral allele.
* `-doSaf 1` and related genotype-likelihood SAF workflows require an ancestral sequence unless simulation input or an alternative major/minor-allele mode is used.
* `-doSaf 5` uses inferred major/minor alleles and therefore does not provide the same polarization as an ancestral-reference analysis.
* Workflows reading precomputed likelihoods or VCF/BCF genotype likelihoods still require compatible allele labels and chromosome coordinates; supplying a file does not remove those assumptions.

Record the ANGSD command, input mode, reference and ancestral files, ploidy, and major/minor-allele strategy with every SAF or SFS result. These choices affect whether the resulting spectrum is folded, polarized, or defined relative to inferred alleles.

Program has a paper
=====
http://www.biomedcentral.com/1471-2105/15/356/abstract
