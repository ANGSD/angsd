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
