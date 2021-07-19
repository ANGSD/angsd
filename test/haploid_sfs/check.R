##this is code for validating that the old -fold 1 in angsd and -dothetas in angsd correspond to the similar analyses that has been moved to realSFS where it is easier to generalize
##this is not meant being run, but will stay in angsd for the time being as documentation

##compare SFS for the folded (in angsd) and fold in realSFS
a <-  scan("output//fold.saf.em.ml")
b <-  scan("output//fold2.saf.em.ml")
b[b!=0]-a
# [1]  2.410800e-02 -2.538800e-02  1.336000e-03 -5.900000e-05  2.000001e-06
# [6]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#[11]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#[16]  0.000000e+00 -2.000000e-06  3.000000e-06 -1.000000e-06 -1.000000e-06
#[21]  9.999999e-07

##compare unfolded popstats
a<-read.table("output/norm.thetas.idx.pestPG")
b<-read.table("output/norm.thetaFromSaf.thetas.idx.pestPG")
a[1]
##                              V1
##1 (0,264255)(1,264256)(0,264256)
b[1]
##                              V1
##1 (0,264255)(1,264256)(0,264256)
as.numeric(a[-1])-as.numeric(b[-1])
## [1]  0.000000e+00  0.000000e+00  0.000000e+00  3.000001e-06  0.000000e+00
## [6] -6.999995e-06 -1.200000e-05  0.000000e+00  0.000000e+00  0.000000e+00
##[11]  0.000000e+00  0.000000e+00  0.000000e+00


a<-read.table("output/norm.thetas.idx.win.pestPG",as.is=T)
b<-read.table("output/norm.thetaFromSaf.thetas.idx.win.pestPG",as.is=T)
table(a[,1]!=b[,1])
##
##FALSE 
##  259 
range(a[,-1]-b[,-1])
##[1] -3e-06  3e-06

##check persite thetas
##../../misc/thetaStat print output/norm.thetas.idx >output/norm.thetas.idx.txt
##../../misc/thetaStat print output/norm.thetaFromSaf.thetas.idx >output/norm.thetaFromSaf.thetas.idx.txt
a<-read.table("output/norm.thetas.idx.txt")
b<-read.table("output/norm.thetaFromSaf.thetas.idx.txt")


#> table(is.infinite(as.matrix(a))==is.infinite(as.matrix(b)))
##
##
##   TRUE 
##1849792 
a<- as.matrix(a)
b<- as.matrix(b)
a[is.infinite(a)] <- 0
b[is.infinite(b)] <- 0
                                        ##
range(a-b)
##[1] -1e-06  1e-06

##Compare unfolded popstats
a<-read.table("output/fold.thetas.idx.pestPG")
b<-read.table("output/fold.thetaFromSaf.thetas.idx.pestPG")
a[1]
##                              V1
##1 (0,264255)(1,264256)(0,264256)
b[1]
##                              V1
##1 (0,264255)(1,264256)(0,264256)

as.numeric(a[-1])-as.numeric(b[-1])
## [1]  0.000000e+00  0.000000e+00  0.000000e+00  3.000001e-06  0.000000e+00
## [6] -6.999995e-06 -1.200000e-05  0.000000e+00  0.000000e+00  0.000000e+00
##[11]  0.000000e+00  0.000000e+00  0.000000e+00


a<-read.table("output/norm.thetas.idx.win.pestPG",as.is=T)
b<-read.table("output/norm.thetaFromSaf.thetas.idx.win.pestPG",as.is=T)
table(a[,1]!=b[,1])
##
##FALSE 
##  259 
range(a[,-1]-b[,-1])
#[1] -3e-06  3e-06

##check persite thetas norm
../../misc/thetaStat print output/norm.thetas.idx >output/norm.thetas.idx.txt
../../misc/thetaStat print output/norm.thetaFromSaf.thetas.idx >output/norm.thetaFromSaf.thetas.idx.txt
../../misc/thetaStat print output/fold.thetas.idx >output/fold.thetas.idx.txt
../../misc/thetaStat print output/fold.thetaFromSaf.thetas.idx >output/fold.thetaFromSaf.thetas.idx.txt


a<-as.matrix(read.table("output/norm.thetas.idx.txt"))
b<-as.matrix(read.table("output/norm.thetaFromSaf.thetas.idx.txt"))
table(a[is.infinite(a)]!=b[is.infinite(b)])

##
## FALSE 
##106721 
a[is.infinite(a)]<-0
b[is.infinite(b)]<-0
 range(a-b)
##[1] -1e-06  1e-06

a<-as.matrix(read.table("output/fold.thetas.idx.txt"))
b<-as.matrix(read.table("output/fold.thetaFromSaf.thetas.idx.txt"))
table(a[is.infinite(a)]!=b[is.infinite(b)])
a[is.infinite(a)]<-0
b[is.infinite(b)]<-0
##
table(a[is.infinite(a)]!=b[is.infinite(b)])

## FALSE 
##792768 

a[is.infinite(a)]<-0
b[is.infinite(b)]<-0
 range(a-b)
##[1] -1e-06  7.2e-05

##per site looks fine, lets look at the global and the window statistic

range(read.table("output/norm.thetas.idx.pestPG")[1,-1]-read.table("output/norm.thetaFromSaf.thetas.idx.pestPG")[1,-1])
##[1] -1.200000e-05  3.000001e-06

range(read.table("output/norm.thetas.idx.win.pestPG")[,-1]-read.table("output/norm.thetaFromSaf.thetas.idx.win.pestPG")[,-1])
#[1] -3e-06  3e-06


range(read.table("output/fold.thetas.idx.pestPG")[1,-1]-read.table("output/fold.thetaFromSaf.thetas.idx.pestPG")[1,-1])
##[1] 0.000000 0.005551
read.table("output/fold.thetas.idx.pestPG")[1,-1]-read.table("output/fold.thetaFromSaf.thetas.idx.pestPG")[1,-1]
##  V2 V3       V4       V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
##1  0  0 0.005551 0.001118  0  0  0  0   0   0   0   0   0

##Difference in 0.005 in watterson and 0.001 in pairwise, that seems reasonable given that the old analyses was with double and the new with floats. This should just be a difference in precision
range(read.table("output/fold.thetas.idx.win.pestPG")[,-1]-read.table("output/fold.thetaFromSaf.thetas.idx.win.pestPG")[,-1])
##[1] -0.000001  0.000181

#again this looks like the same

