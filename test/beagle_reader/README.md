# changelog beagle reading [beagleReader.cpp]
## 220528 isinaltinkaya

Changed beagle field parsing from strtok_r to strrchr to handle cases where multiple _s are used.

```
../../../fixed_angsd/angsd/angsd -doMaf 4 -beagle test1.beagle -fai ref1.fa.fai -out test1
../../../fixed_angsd/angsd/angsd -doMaf 4 -beagle test2.beagle -fai ref2.fa.fai -out test2
../../../fixed_angsd/angsd/angsd -doMaf 4 -beagle test3.beagle -fai ref3.fa.fai -out test3
/science/willerslev/scratch/pfs488/angsd/issue388/main_angsd/angsd/angsd -doMaf 4 -beagle test1.beagle -fai ref1.fa.fai -out main1
/science/willerslev/scratch/pfs488/angsd/issue388/main_angsd/angsd/angsd -doMaf 4 -beagle test2.beagle -fai ref2.fa.fai -out main2
```

```
==> test1.beagle <==
1_5 c a 0.9 0.05 0.05 0.8 0.10 0.10 0.9 0.05 0.05
1_6 g t 0.9 0.05 0.05 0.8 0.10 0.10 0.9 0.05 0.05
1_7 g c 0.9 0.05 0.05 0.6 0.20 0.20 0.8 0.10 0.10
2_10 t a 0.9 0.05 0.05 0.9 0.05 0.05 0.9 0.05 0.05
2_11 c a 0.9 0.05 0.05 0.9 0.05 0.05 0.9 0.05 0.05
2_12 g c 0.9 0.05 0.05 0.9 0.05 0.05 0.7 0.15 0.15

==> test2.beagle <==
chr1_5 c a 0.9 0.05 0.05 0.8 0.10 0.10 0.9 0.05 0.05
chr1_6 g t 0.9 0.05 0.05 0.8 0.10 0.10 0.9 0.05 0.05
chr1_7 g c 0.9 0.05 0.05 0.6 0.20 0.20 0.8 0.10 0.10
chr2_10 t a 0.9 0.05 0.05 0.9 0.05 0.05 0.9 0.05 0.05
chr2_11 c a 0.9 0.05 0.05 0.9 0.05 0.05 0.9 0.05 0.05
chr2_12 g c 0.9 0.05 0.05 0.9 0.05 0.05 0.7 0.15 0.15

==> test3.beagle <==
chr_1_5 c a 0.9 0.05 0.05 0.8 0.10 0.10 0.9 0.05 0.05
chr_1_6 g t 0.9 0.05 0.05 0.8 0.10 0.10 0.9 0.05 0.05
chr_1_7 g c 0.9 0.05 0.05 0.6 0.20 0.20 0.8 0.10 0.10
chr_2_10 t a 0.9 0.05 0.05 0.9 0.05 0.05 0.9 0.05 0.05
chr_2_11 c a 0.9 0.05 0.05 0.9 0.05 0.05 0.9 0.05 0.05
chr_2_12 g c 0.9 0.05 0.05 0.9 0.05 0.05 0.7 0.15 0.15
```

```
==> ref1.fa.fai <==
1       249250621       3       50      51

==> ref2.fa.fai <==
chr1    249250621       3       50      51

==> ref3.fa.fai <==
chr_1   249250621       3       50      51
```


```
$ zcat *mafs.gz
chromo  position        major   minor   PPmaf   nInd
1       6       G       T       0.100000        3
1       7       G       C       0.175000        3
2       10      T       A       0.075000        3
2       11      C       A       0.075000        3
2       12      G       C       0.125000        3
chromo  position        major   minor   PPmaf   nInd
chr1    6       G       T       0.100000        3
chr1    7       G       C       0.175000        3
chr2    10      T       A       0.075000        3
chr2    11      C       A       0.075000        3
chr2    12      G       C       0.125000        3
chromo  position        major   minor   PPmaf   nInd
1       6       G       T       0.100000        3
1       7       G       C       0.175000        3
2       10      T       A       0.075000        3
2       11      C       A       0.075000        3
2       12      G       C       0.125000        3
chromo  position        major   minor   PPmaf   nInd
chr1    6       G       T       0.100000        3
chr1    7       G       C       0.175000        3
chr2    10      T       A       0.075000        3
chr2    11      C       A       0.075000        3
chr2    12      G       C       0.125000        3
chromo  position        major   minor   PPmaf   nInd
chr_1   6       G       T       0.100000        3
chr_1   7       G       C       0.175000        3
chr_2   10      T       A       0.075000        3
chr_2   11      C       A       0.075000        3
chr_2   12      G       C       0.125000        3
```

### Note:
[BUG?] main angsd and current version skips line1.

[TODO] add this change to tests after possible bug above is resolved.
