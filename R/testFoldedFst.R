library(SQUAREM)
library(parallel)

norm <- function(x,...) x/sum(x,...)
na2zero <-function(x,...) {x[is.na(x)]<-0;return(x)}
remapper <- function(x,y){
    m<-matrix(1:(x*y),x)
    m1<-m
    nal=x+y-2
    for( i in 1:x)
        for(j in 1:y){
            af<-(i+j-2)/nal
            if(af>0.5)
                m1[i,j]<-m[x-i+1,y-j+1]
        }
    colnames(m1)<-0:(y-1)
    rownames(m1)<-0:(x-1)
    m1
}




##This code is modified from ANGSD R script, which estimates Reynolds 1983 Fst, to estimate Hudsons 1992 Fst as interpreted by Bhatia 2013 ###
getFst<-function(est){
    est[is.na(est)] <- 0
    N1<-nrow(est)-1
    N2<-ncol(est)-1
    ##    cat("N1: ",N1 ," N2: ",N2,"\n")
    est0<-est
    est0[1,1]<-0
    est0[N1+1,N2+1]<-0
    est0<-est0/sum(est0)

    aMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    aMat.ss<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    baMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    for(a1 in 0:(N1))
        for(a2 in 0:(N2)){
            p1 <- a1/N1
            p2 <- a2/N2
            q1 <- 1 - p1
            q2 <- 1 - p2
            N <- (p1-p2)^2
            D <- p1*(1-p2)+p2*(1-p1)
            aMat[a1+1,a2+1]<<-N
            baMat[a1+1,a2+1]<<-D
            #sample size correction
            N.ss <- (p1-p2)^2-((p1*(1-p1))/(N1-1))-((p2*(1-p2))/(N2-1))
            aMat.ss[a1+1,a2+1]<<-N.ss
        }
    ## sample size corrected moment estimator
    ss <- sum(est0*aMat.ss,na.rm=T)/sum(est0*baMat,na.rm=T) 
    c(fstSS=ss)
}

#######################HERE IS THE ORIGINAL (ANGSD) REYNOLD's ESTIMATOR#########################################

getReynoldsFst<-function(est){
    est[is.na(est)] <- 0
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

fold2d <- function(x){
    thMid <- (nrow(x)+ncol(x)-2)/2
    ##   ret <- matrix(NA,min(thMid+1,nrow(x)),min(thMid+1,ncol(x)))
    ret <- matrix(NA,nrow(x),ncol(x))
    colnames(ret) <- 0:(ncol(ret)-1)
    rownames(ret) <- 0:(nrow(ret)-1)
    for(i in 1:nrow(x))
        for(j in 1:ncol(x)){
            af<- ((i-1+j-1)/sum(dim(x)-1))
      ##      cat("i:",i,"j:",j,"af:",af);
            if( ((i-1+j-1)/sum(dim(x)-1))<=0.5){
        ##        cat('\n')
                ret[i,j] <- ifelse(is.na(ret[i,j]),x[i,j],ret[i,j]+x[i,j])
            }else{
          ##      cat("AF:",sum(dim(x)-1),"af:",af,"(",i,",",j,")(",nrow(x)-i+1,",",ncol(x)-j+1,")\n")
                ret[nrow(x)-i+1,ncol(x)-j+1] <- ifelse(is.na(ret[nrow(x)-i+1,ncol(x)-j+1]),x[i,j],ret[nrow(x)-i+1,ncol(x)-j+1]+x[i,j])
            }
        }
    if(abs(sum(x)-sum(ret,na.rm=T))>1e-12)
        stop("Problem with difference of values:",sum(ret)," vs ",sum(x))
    return(ret)
}


if(FALSE){
    ##data simulated like
    ###MS=/willerslev/software/msms/bin/msms
    ##$MS -ms 90 10 -t 930 -r 400 -I 3 26 30 34 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1 |bgzip -c >${MSMS}
    ## see angsd/test/testFst_folded.sh for commandline stuff
    
    source('testFoldedFst.R')
    source("readms.output.R")
    nPop1 <- 26
    nPop2 <- 30 
    nPop3 <- 34
    a <- read.ms.output(file="../test/fst_folded/input/msout")
    p1.d <- unlist((sapply(a$gam,function(x) colSums(x[1:nPop1,]))))
    p2.d <- unlist((sapply(a$gam,function(x) colSums(x[(nPop1+1):(nPop1+nPop2),]))))
    p3.d <- unlist((sapply(a$gam,function(x) colSums(x[(nPop1+nPop2+1):(nPop1+nPop2+nPop3),]))))
    p1.p2.sfs <-  t(sapply(0:nPop1,function(x) table(factor(p2.d[p1.d==x],levels=0:nPop2))))
    p1.p3.sfs <-  t(sapply(0:nPop1,function(x) table(factor(p3.d[p1.d==x],levels=0:nPop3))))
    p2.p3.sfs <-  t(sapply(0:nPop2,function(x) table(factor(p3.d[p2.d==x],levels=0:nPop3))))
    a.p1.sfs <- scan("../test/fst_folded/output/pop1.saf.idx.ml")
    a.p2.sfs <- scan("../test/fst_folded/output/pop2.saf.idx.ml")
    a.p3.sfs <- scan("../test/fst_folded/output/pop3.saf.idx.ml")


    par(mfrow=c(3,1))
    barplot(rbind(table(p1.d),a.p1.sfs),be=T)
    barplot(rbind(table(p2.d),a.p2.sfs),be=T)
    barplot(rbind(table(p3.d),a.p3.sfs),be=T)
    ##All good sofar with the 1dsfs

    a.p1.p2.sfs <- matrix(scan("../test/fst_folded/output/pop1.pop2.saf.idx.ml"),byrow=T,nrow(p1.p2.sfs),ncol(p1.p2.sfs))
    a.p1.p3.sfs <- matrix(scan("../test/fst_folded/output/pop1.pop3.saf.idx.ml"),byrow=T,nrow(p1.p3.sfs),ncol(p1.p3.sfs))
    a.p2.p3.sfs <- matrix(scan("../test/fst_folded/output/pop2.pop3.saf.idx.ml"),byrow=T,nrow(p2.p3.sfs),ncol(p2.p3.sfs))
    par(mfrow=c(3,1))
    hist(p1.p2.sfs-a.p1.p2.sfs)
    hist(p1.p3.sfs-a.p1.p3.sfs)
    hist(p2.p3.sfs-a.p2.p3.sfs)
    par(mfrow=c(3,1))
    plot(p1.p2.sfs,a.p1.p2.sfs)
    plot(p1.p3.sfs,a.p1.p3.sfs)
    plot(p2.p3.sfs,a.p2.p3.sfs)

    ##All good sofar with the 2dsfs. , we are only simulatre 64k sites with 837 categories of the spectra

    ##now lets redo with folding
    a.p1.p2.sfs.fold <- matrix(scan("../test/fst_folded/output/pop1.pop2.saf.idx.ml.fold"),byrow=T,nrow(p1.p2.sfs),ncol(p1.p2.sfs))
    a.p1.p3.sfs.fold <- matrix(scan("../test/fst_folded/output/pop1.pop3.saf.idx.ml.fold"),byrow=T,nrow(p1.p3.sfs),ncol(p1.p3.sfs))
    a.p2.p3.sfs.fold <- matrix(scan("../test/fst_folded/output/pop2.pop3.saf.idx.ml.fold"),byrow=T,nrow(p2.p3.sfs),ncol(p2.p3.sfs))
    par(mfrow=c(3,1))
    hist(na2zero(fold2d(p1.p2.sfs))-a.p1.p2.sfs.fold)
    hist(na2zero(fold2d(p1.p3.sfs))-a.p1.p3.sfs.fold)
    hist(na2zero(fold2d(p2.p3.sfs))-a.p2.p3.sfs.fold)
    par(mfrow=c(3,1))
    plot(na2zero(fold2d(p1.p2.sfs)),a.p1.p2.sfs.fold)
    plot(na2zero(fold2d(p1.p3.sfs)),a.p1.p3.sfs.fold)
    plot(na2zero(fold2d(p2.p3.sfs)),a.p2.p3.sfs.fold)
    par(mfrow=c(3,1))
    boxplot(na2zero(fold2d(p1.p2.sfs))-a.p1.p2.sfs.fold)
    boxplot(na2zero(fold2d(p1.p3.sfs))-a.p1.p3.sfs.fold)
    boxplot(na2zero(fold2d(p2.p3.sfs))-a.p2.p3.sfs.fold)

    
    ##This looks fine, but lets implement a glf reader in R and implement the folded optimization so we are sure things really work
    
    ##code below implements the stuff that happens in realSFS.cpp
    readdata<-function(fname,x,nsites){
        ff <- gzfile(fname,"rb")
        seek(ff,8)
        m<-matrix(readBin(ff,double(),x*nsites,size=4),ncol=x,byrow=TRUE)
        close(ff)
        return(m)
    }
    gl1<-t(exp(readdata('../test/fst_folded/output/pop1.saf.gz',nPop1+1,1e9)))
    gl2<-t(exp(readdata('../test/fst_folded/output/pop2.saf.gz',nPop2+1,1e9)))
    gl12<-matrix(NA,nrow=(nPop1+1)*(nPop2+1),ncol=ncol(gl1))
    for(i in 1:ncol(gl1))
        gl12[,i] <- as.numeric(t(gl1[,i] %*% t(gl2[,i])))
        

    mylik<-function(x,d)
        sum(log(colSums(x*d)))
    emStep <- function(x,d)
        rowMeans(apply(x*d,2,norm))
    ##start is in log ,d is in log
    emer<-function(start,d,niter=100,tole=1e-6){
        if(length(start)!=nrow(d))
            stop("")
        oldlik<-mylik(start,d)
        cat('startlik: ',oldlik,"\n")
        for(i in 1:niter){
            start<-emStep(start,d)
            lik<-sum(log(colSums(start*d)))
            cat("niter: ",i," lik: ",lik, 'dif: ',abs(oldlik-lik),"\n")
            if(abs(oldlik-lik)<tole)
                break
            oldlik<-lik
        }
        return(start)
    }
  
    sfs.2d.est <-squarem(norm(runif( (nPop1+1)*(nPop2+1))),emStep,mylik,d=gl12)
    plot(sfs.2d.est$par*ncol(gl12),scan("../test/fst_folded/output/pop1.pop2.saf.idx.ml"))
    plot(sfs.2d.est$par*ncol(gl12)-scan("../test/fst_folded/output/pop1.pop2.saf.idx.ml"))
    hist(sfs.2d.est$par*ncol(gl12)-scan("../test/fst_folded/output/pop1.pop2.saf.idx.ml"),br=101)
    ##these estimates looks like the c implemenation
    

    gl12.fold <- mclapply(1:ncol(gl1),function(i) as.numeric(t(fold2d(gl1[,i] %*% t(gl2[,i]))) ),mc.cores=60)
    gl12.fold <- matrix(unlist(gl12.fold),nrow=(nPop1+1)*(nPop2+1))
    gl12.fold[is.na(gl12.fold)] <- 0
    sfs.2d.fold.est <-squarem(norm(runif((nPop1+1)*(nPop2+1))),emStep,mylik,d=gl12.fold)
    res.fold<-matrix(sfs.2d.fold.est$par,nPop1+1,nPop2+1,byrow=T)*ncol(gl12.fold)
    tmp<-fold2d(p1.p2.sfs)
    tmp[is.na(tmp)] <- 0
    plot(res.fold,tmp)
    ##seems to work
    
    start<-norm(fold2d(matrix(1:(nPop1+1)*(nPop2+1),nPop1+1,nPop2+1,byrow=T)),na.rm=T)
    start[is.na(start)] <- 0
    start<-as.numeric(t(start))
    sfs.2d.fold.est2 <-squarem(start,emStep,mylik,d=gl12.fold,control=list(maxiter=1500))
    hist(sfs.2d.fold.est2$par*ncol(gl12.fold)-scan('../test/fst_folded/output/pop1.pop2.saf.idx.ml.fold'))
    ##also works now. lets stop here its 2.13am
    plot(sfs.2d.fold.est2$par*ncol(gl12.fold)-scan('../test/fst_folded/output/pop1.pop2.saf.idx.ml.fold'))

    ##lets check the fst values for the unfolded
    
}
