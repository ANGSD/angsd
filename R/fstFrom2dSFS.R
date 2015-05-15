
reynolds<-function(pl1,pl2,n1,n2) { ##from matteo
    somma=sommaden=0;
    alfa1=1-((pl1^2)+((1-pl1)^2))
    alfa2=1-((pl2^2)+((1-pl2)^2))
    Al = (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) - (((n1+n2)*(n1*alfa1+n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))
    AlBl= (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) + (((4*n1*n2 - n1 - n2)*(n1*alfa1 + n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))
    if (!is.na(Al) & !is.na(AlBl)) {
        somma=somma+Al
        sommaden=sommaden+AlBl
    }
    if (somma==0 & sommaden==0) {
	reyn=NA
    } else {
	reyn=somma/sommaden
 	if(reyn<0) reyn=0
    }
    reyn
}

getFst<-function(est){
    N1<-nrow(est)-1
    N2<-ncol(est)-1
    est0<-est
    est0[1,1]<-0
    est0[N1+1,N2+1]<-0
    est0<-est0/sum(est0)
    
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
            
            al <-  1/2 * ( (p1-p2)^2 + (q1-q2)^2) -
                (N1+N2)*        (N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
            bal <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) +
                (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
            
            aMat[a1+1,a2+1]<-al
            baMat[a1+1,a2+1]<-bal
            ##  print(signif(c(a1=a1,a2=a2,p1=p1,p2=p2,al1=alpha1,al2=alpha2,al),2))
        }
    # unweighted average of single-locus ratio estimators
 fstU <-   sum(est0*(aMat/baMat),na.rm=T)
    # weighted average of single-locus ratio estimators
 fstW <-   sum(est0*aMat,na.rm=T)/sum(est0*baMat,na.rm=T)
c(fstW=fstW,fstU=fstU)
}
##
##matteo
getFst2 <- function(sfs){
    sfs[1,1]=NA
    sfs[nrow(sfs),ncol(sfs)]=NA
    
    RELATIVE=1
    if (max(sfs, na.rm=T)<=1)
        RELATIVE=0
    if (RELATIVE)
        sfs=sfs/sum(sfs,na.rm=T)
    
    nind1=(nrow(sfs)-1)/2
    nind2=(ncol(sfs)-1)/2
    
    nind1=(nrow(sfs)-1)/2
    nind2=(ncol(sfs)-1)/2
    ##cat("Pop1 (rows) has",nind1,"individuals")
    ##cat("\nPop2 (cols) has",nind2,"individuals")
    
    fsts=matrix(NA, ncol=ncol(sfs), nrow=nrow(sfs))
    for (i in 1:nrow(sfs)) {
        for (j in 1:ncol(sfs)) {
            f1=(i-1)/(nind1*2); f2=(j-1)/(nind2*2)
                                        #cat("\n",f1," ",f2)
            fsts[i,j]=reynolds(f1,f2,nind1*2,nind2*2)  
        }
    }
    
                                        #cat("\nsum sfs:",sum(sfs))
        
    if (RELATIVE) {
        fst=sum(sfs*fsts,na.rm=T)
    } else {
        fst=sum(sfs*fsts,na.rm=T)/sum(sfs)
    }
    
    c(fstM=fst)
}



if(FALSE){
    if(FALSE){
        ##simulate data with msms
       nRep <- 10
       nPop1 <- 24
       nPop2 <- 16
       cmd <- paste("msms -ms",nPop1+nPop2,nRep,"-t 930 -r 400 -I 2",nPop1,nPop2,"0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt ",sep=" ")
       system(cmd)
       ##system("msms -ms 40 1 -t 930 -r 400 -I 2 20 20 0 -g 1 9.70406 -n 1 2 -n 2 1 -ma x 0.0 0.0 x -ej 0.07142857 2 1  >msoutput.txt  ")
       
       source("../R/readms.output.R")
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
        source("fstFrom2dSFS.R")
        ##type1
        getFst(sfs.2d)
        ##type2
        getFst2(sfs.2d)
    }
}
