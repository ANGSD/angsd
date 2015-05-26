getFst<-function(est){
    N1<-nrow(est)-1
    N2<-ncol(est)-1
    cat("N1: ",N1 ," N2: ",N2,"\n")
    est0<-est
    est0[1,1]<-0
    est0[N1+1,N2+1]<-0
    est0<-est0/sum(est0)
    
    aMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    baMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    for(a1 in 0:(N1))
        for(a2 in 0:(N2)){
            p1 <- a1/N1
            p2 <- a2/N2
            q1 <- 1 - p1
            q2 <- 1 - p2
            alpha1 <- 1 - (p1^2 + q1^2)
            alpha2 <- 1 - (p2^2 + q2^2)
            
            al <-  1/2 * ( (p1-p2)^2 + (q1-q2)^2) - (N1+N2) *  (N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
            bal <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) + (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
            aMat[a1+1,a2+1]<<-al
            baMat[a1+1,a2+1]<<-bal
            ##  print(signif(c(a1=a1,a2=a2,p1=p1,p2=p2,al1=alpha1,al2=alpha2,al),2))
        }
    ## unweighted average of single-locus ratio estimators
    fstU <-   sum(est0*(aMat/baMat),na.rm=T)
    ## weighted average of single-locus ratio estimators
    fstW <-   sum(est0*aMat,na.rm=T)/sum(est0*baMat,na.rm=T)
    c(fstW=fstW,fstU=fstU)
}

getCoefs <- function(N1,N2){
    cat("N1: ",N1 ," N2: ",N2,"\n")
    aMat<-matrix(NA,nrow=N1+1,ncol=N2+1)
    baMat<-matrix(NA,nrow=N1+1,ncol=N2+1)
    for(a1 in 0:(N1))
        for(a2 in 0:(N2)){
            p1 <- a1/N1
            p2 <- a2/N2
            q1 <- 1 - p1
            q2 <- 1 - p2
            alpha1 <- 1 - (p1^2 + q1^2)
            alpha2 <- 1 - (p2^2 + q2^2)
            
            al <-  1/2 * ( (p1-p2)^2 + (q1-q2)^2) - (N1+N2) *  (N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
            bal <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) + (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
            aMat[a1+1,a2+1]<-al
            baMat[a1+1,a2+1]<-bal
           
        }
    list(a1=aMat,b1=baMat)
}

getFst(sfs.2d)



if(FALSE){
    if(FALSE){
        ##simulate data with msms
       nRep <- 100
       nPop1 <- 24
       nPop2 <- 16
       cmd <- paste("msms -ms",nPop1+nPop2,nRep,"-t 930 -r 400 -I 2",nPop1,nPop2,"0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt ",sep=" ")
       system(cmd)
       ##system("msms -ms 40 1 -t 930 -r 400 -I 2 20 20 0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt  ")
       
       source("../R/readms.output.R")
   }
    to2dSFS <- function(p1.d,p2.d,nPop1,nPop2)
        sapply(0:nPop1,function(x) table(factor(p2.d[p1.d==x],levels=0:nPop2)))
   

    if(FALSE){
        ##use R to calculate SFS for each pop and 2dsfs
        source("../R/readms.output.R")
        a<- read.ms.output(file="msoutput.txt")
        
        p1.d <- unlist((sapply(a$gam,function(x) colSums(x[1:nPop1,]))))
        p2.d <- unlist((sapply(a$gam,function(x) colSums(x[-c(1:nPop1),]))))
        par(mfrow=c(1,2))
        barplot(table(p1.d))
        barplot(table(p2.d))
        sfs.2d <- t(sapply(0:nPop1,function(x) table(factor(p2.d[p1.d==x],levels=0:nPop2))))
        sfs.2d.sub1 <- to2dSFS(p1.d[c(1:40e3)],p2.d[1:40e3],nPop1,nPop2)
        sfs.2d.sub2 <- to2dSFS(p1.d[-c(1:40e3)],p2.d[-c(1:40e3)],nPop1,nPop2)
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
       system("../misc/realSFS pop1.saf.idx pop2.saf.idx -maxIter 500 -p 20  >pop1.pop2.saf.idx.ml")
   }
            
    
    if(FALSE){
        source("fstFrom2dSFS.R")
        ##type1
        getFst(sfs.2d)
    }
    if(FALSE){
        norm <- function(x) x/sum(x)
        prior <- norm(sfs.2d)
        coef <- getCoefs(nrow(prior)-1,ncol(prior)-1)
        ab <- matrix(NA,ncol=2,nrow=length(p1.d))
        for(r in 1:length(p1.d)){
            if((r %% 100 )==0 )
                cat( "\r  ",r,"/ ", length(p1.d))
            mat <- matrix(0,ncol=nPop1+1,nrow=nPop2+1)
            mat[p2.d[r]+1,p1.d[r]+1] <- 1
        #    mat <- mat*prior
            ab[r,] <- c(sum(mat*coef$a1),sum(mat*coef$b1))
        }
        
    }
    if(FALSE){
        est <- matrix(as.integer(scan("pop1.pop2.saf.idx.ml")),byrow=T,ncol=nPop2+1)
        

    }
    
    
}
