## [![Build Status](https://github.com/ANGSD/angsd/actions/workflows/build-tests.yml/badge.svg)](https://github.com/ANGSD/angsd/actions/workflows/build-tests.yml) 


Program for analysing NGS data.
=====
 

Documentation: https://github.com/ANGSD/angsd/wiki
 
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

3) Let make fetch htslib (simplest)

```
git clone https://github.com/angsd/angsd.git;
cd angsd; make
```
`make` clones the newest htslib `develop` into `angsd/htslib` and updates it on later runs (offline it keeps the copy it has). Release tarballs ship with htslib included, so `make` builds them without network access.

Versions and releases
=====
* Releases are annotated tags on master, e.g. `0.942`, and each one gets a GitHub release with a source tarball (made by `.github/workflows/release.yml`).
* `angsd` reports its version from `git describe`: `0.942` on a release, `0.942-3-gabc1234` three commits after it, with `-dirty` if the source has uncommitted changes. Source tarballs carry the version in the `VERSION` file.
* To release: `git tag -a 0.943 -m "angsd 0.943" && git push origin 0.943`.



Notes
====
* I've switched over to using htslib for parsing single reads (to allow for CRAM reading, while avoid having to write my own CRAM parser). I'm still using my own readpools. Users should therefore also download and install htslib.

Reproducibility record
=====
For each analysis, retain the following alongside the output files:

* the complete command line and the generated `.args` file;
* the ANGSD version, build date and HTSlib version;
* the BAM/VCF/BCF input list and the reference and ancestral sequence identifiers;
* ploidy, genotype-likelihood and posterior-probability modes;
* site, depth, mapping-quality, base-quality and SNP filters;
* the major/minor-allele and folded/unfolded SFS choices;
* the random seed, thread count and any region or site-index files.

ANGSD defaults and output conventions can change between releases. A results file without this information is difficult to reproduce, even when the same nominal command is reused later.

Program has a paper
=====
http://www.biomedcentral.com/1471-2105/15/356/abstract
