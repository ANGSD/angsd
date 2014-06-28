if(FALSE){
  ##generate data DONT RUN
  nRep <- 10
  nPop1 <- 24
  nPop2 <- 16
  cmd <- paste("msms -ms",nPop1+nPop2,nRep,"-t 930 -r 400 -I 2",nPop1,nPop2,"0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt ",sep=" ")
  system(cmd)
  ##system("msms -ms 40 1 -t 930 -r 400 -I 2 20 20 0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt  ")
  
  source("readms.output.R")

  a<- read.ms.output(file="msoutput.txt")

  p1.d <- unlist((sapply(a$gam,function(x) colSums(x[1:nPop1,]))))
  p2.d <- unlist((sapply(a$gam,function(x) colSums(x[-c(1:nPop1),]))))
  par(mfrow=c(1,2))
  barplot(table(p1.d))
  barplot(table(p2.d))
  
  sfs.2d <- sapply(0:nPop1,function(x) table(factor(p2.d[p1.d==x],levels=0:nPop2)))

  ##begin cmdline
  system("./msToGlf6 -in msoutput.txt -out raw -singleOut 1 -regLen 0 -depth 8 -err 0.005")
  system("./angsd0.540/angsd -sim1 raw.glf.gz -nInd 20 -realSFS 1 -out rawa -from 0 -to 11 -nThreads 20")
  system("./angsd0.540/angsd -sim1 raw.glf.gz -nInd 20 -realSFS 1 -out rawb -from 12 -to 19 -nThreads 20")
  system("./emOptim2 2dsfs rawb 24 16 -tole 1e-6  -P 10 >resVars")


  tsk1 <- norm(exp(t(read.table("resVars"))))
  range(norm(sfs.2d)-tsk1)
  ##[1] -0.0005869855  0.0005790802

  barplot(rbind(rowSums(norm(sfs.2d)),rowSums(tsk1)))
  barplot(rbind(colSums(norm(sfs.2d)),colSums(tsk1)))
  

  system("./msToGlf6 -in msoutput.txt -out rawInvar -singleOut 1 -regLen 1000000 -depth 8 -err 0.005")
  system("./angsd -sim1 rawInvar.glf.gz -nInd 20 -realSFS 1 -out rawaInvar -from 0 -to 11 -nThreads 20")
  system("./angsd -sim1 rawInvar.glf.gz -nInd 20 -realSFS 1 -out rawbInvar -from 12 -to 19 -nThreads 20")
  system("./emOptim2 2dsfs rawbInvar 24 16 -tole 1e-6  -P 24 >resInVars")##tolerence achieved after 180 iterations

  tsk2 <- norm(exp(t(read.table("resInVars"))))
  range(norm(sfs.2d[-1,-1])-norm(tsk2[-1,-1]))
  ##[1] -0.0008945658  0.0106200507

  barplot(rbind(rowSums(norm(sfs.2d[-1,-1])),rowSums(norm(tsk2[-1,-1]))))
  barplot(rbind(colSums(norm(sfs.2d[-1,-1])),colSums(norm(tsk2[-1,-1]))))

}

