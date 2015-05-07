##angsd new: 	->-> angsd version: 0.801-51-g156039a (htslib: 1.2.1-69-gb79f40a) build(May  7 2015 15:31:53)
##angsd old: 	-> angsd version: 0.801-27-ga699b44 (htslib: 1.2.1-69-gb79f40a) build(May  7 2015 15:30:06)

norm <- function(x) x/sum(x)
if(FALSE){
  ##generate data DONT RUN
   if(FALSE){
       ##simulate data with msms
       nRep <- 10
       nPop1 <- 24
       nPop2 <- 16
       cmd <- paste("msms -ms",nPop1+nPop2,nRep,"-t 930 -r 400 -I 2",nPop1,nPop2,"0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt ",sep=" ")
       system(cmd)
       ##system("msms -ms 40 1 -t 930 -r 400 -I 2 20 20 0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt  ")
       
       source("readms.output.R")
   }
   if(FALSE){
       ##use R to calculate SFS for each pop and 2dsfs
       source("../R/readms.output.R")
       a<- read.ms.output(file="msoutput.txt")
       
       p1.d <- unlist((sapply(a$gam,function(x) colSums(x[1:nPop1,]))))
       p2.d <- unlist((sapply(a$gam,function(x) colSums(x[-c(1:nPop1),]))))
       par(mfrow=c(1,2))
       barplot(table(p1.d))
       barplot(table(p2.d))
       
       sfs.2d <- sapply(0:nPop1,function(x) table(factor(p2.d[p1.d==x],levels=0:nPop2)))
   }
   if(FALSE){
       ##generate ANGSD inputfiles without invariable sites and run it
       system("../misc/msToGlf -in msoutput.txt -out raw -singleOut 1 -regLen 0 -depth 8 -err 0.005")
       system("../misc/splitgl raw.glf.gz 20 1 12 >pop1.glf.gz")
       system("../misc/splitgl raw.glf.gz 20 13 20 >pop2.glf.gz")
       system("echo \"1 250000000\" >fai.fai")
       system("../angsd -glf pop1.glf.gz -nind 12 -doSaf 1 -out pop1 -fai fai.fai -issim 1")
       system("../angsd -glf pop2.glf.gz -nind 8 -doSaf 1 -out pop2 -fai fai.fai -issim 1")
       system("../misc/realSFS pop1.saf.idx >pop1.saf.idx.ml")
       system("../misc/realSFS pop2.saf.idx >pop2.saf.idx.ml")
       system("../misc/realSFS pop1.saf.idx pop2.saf.idx >pop1.pop2.saf.idx.ml")
   }
   if(FALSE){
       pop1 <- exp(scan("pop1.saf.idx.ml"))
       pop2 <- exp(scan("pop2.saf.idx.ml"))
       pop1.pop2 <- matrix(exp(scan("pop1.pop2.saf.idx.ml")),nPop1+1,byrow=T)
       par(mfrow=c(1,2))
       barplot(rbind(norm(table(p1.d)),pop1),be=T,main="only varsites pop1")
       barplot(rbind(norm(table(p2.d)),pop2),be=T,main="only varsites pop2")
       range(norm(sfs.2d)-t(pop1.pop2))
       ##[1] -0.0005150658  0.0004685074

       barplot(rbind(rowSums(norm(t(sfs.2d))),rowSums(pop1.pop2)),be=T)
       barplot(rbind(colSums(norm(t(sfs.2d))),colSums(pop1.pop2)),be=T)

   }
   if(FALSE){
       ##simulate angsd inputfiles with invariable sites and run it
       system("../misc/msToGlf -in msoutput.txt -out raw -singleOut 1 -regLen 10000000 -depth 8 -err 0.005")
       system("../misc/splitgl raw.glf.gz 20 1 12 >pop1.glf.gz")
       system("../misc/splitgl raw.glf.gz 20 13 20 >pop2.glf.gz")
       system("echo \"1 250000000\" >fai.fai")
       system("../angsd -glf pop1.glf.gz -nind 12 -doSaf 1 -out pop1 -fai fai.fai -issim 1")
       system("../angsd -glf pop2.glf.gz -nind 8 -doSaf 1 -out pop2 -fai fai.fai -issim 1")
       system("../misc/realSFS pop1.saf.idx >pop1.saf.idx.ml")
       system("../misc/realSFS pop2.saf.idx >pop2.saf.idx.ml")
       system("../misc/realSFS pop1.saf.idx pop2.saf.idx >pop1.pop2.saf.idx.ml")

   }
   if(FALSE){
       pop1 <- norm(exp(scan("pop1.saf.idx.ml"))[-1])
       pop2 <- norm(exp(scan("pop2.saf.idx.ml"))[-1])
       pop1.pop2 <- matrix(exp(scan("pop1.pop2.saf.idx.ml")),nPop1+1,byrow=T)
       par(mfrow=c(1,2))
       barplot(rbind(norm(table(p1.d)[-1]),pop1),be=T,main="varsites pop1")
       barplot(rbind(norm(table(p2.d)[-1]),pop2),be=T,main="varsites pop2")
       pop1.pop2[1,1] <- 0
       pop1.pop2[nrow(pop1.pop2),ncol(pop1.pop2)] <- 0
       pop1.pop2 <- norm(pop1.pop2)
       range(norm(sfs.2d)-t(pop1.pop2))
#[1] -0.0005047007  0.0007018686
 
       barplot(rbind(rowSums(norm(t(sfs.2d))),rowSums(pop1.pop2)),be=T)
       barplot(rbind(colSums(norm(t(sfs.2d))),colSums(pop1.pop2)),be=T)

   }
   if(FALSE){
       ##just redo the angsd and optimization
       ##git checkout master ;make clean;make
       system("../angsd -glf pop1.glf.gz -nind 12 -doSaf 1 -out pop1 -fai fai.fai -issim 1 -P 10")
       system("../angsd -glf pop2.glf.gz -nind 8 -doSaf 1 -out pop2 -fai fai.fai -issim 1 -P 10")
       system("../misc/realSFS pop1.saf 24 -P 60 >pop1.saf.idx.ml")
       system("../misc/realSFS pop2.saf 16 -P 60 >pop2.saf.idx.ml")
       system("../misc/realSFS 2dsfs pop1.saf pop2.saf 24 16 -P 60 >pop1.pop2.saf.idx.ml")
   }
   if(FALSE){
       pop1 <- norm(exp(scan("pop1.saf.idx.ml"))[-1])
       pop2 <- norm(exp(scan("pop2.saf.idx.ml"))[-1])
       pop1.pop2 <- matrix(exp(scan("pop1.pop2.saf.idx.ml")),nPop1+1,byrow=T)
       par(mfrow=c(1,2))
       barplot(rbind(norm(table(p1.d)[-1]),pop1),be=T,main="varsites pop1")
       barplot(rbind(norm(table(p2.d)[-1]),pop2),be=T,main="varsites pop2")
       pop1.pop2[1,1] <- 0
       pop1.pop2[nrow(pop1.pop2),ncol(pop1.pop2)] <- 0
       pop1.pop2 <- norm(pop1.pop2)
       range(norm(sfs.2d)-t(pop1.pop2))
       ##       [1] -0.0005046856  0.0007019217

   }
   
}

