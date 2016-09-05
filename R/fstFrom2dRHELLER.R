##This code is modified from ANGSD R script, which estimates Reynolds 1983 Fst, to estimate Hudsons 1992 Fst as interpreted by Bhatia 2013 ###

##read in full 2dsfs:
full<-scan("2dsfs_auto_ancCorrect_indFilters_chapmani_MB.sfs")
fullM<-matrix(full, ncol=35, nrow=17, byrow=T)
##read in subsampled 2dsfs:
r3<-scan("2dsfs_ancCorrect_random3_chapmani_MB.sfs")
m3<-matrix(r3, ncol=7, nrow=7, byrow=T)

##define Fst function
getFst<-function(est){
    N1<-nrow(est)-1
    N2<-ncol(est)-1
    cat("N1: ",N1 ," N2: ",N2,"\n")
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
