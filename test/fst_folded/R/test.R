norm <- function(x,...) x/sum(x,...)
na2zero <-function(x,...) {x[is.na(x)]<-0;return(x)}

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



npop1 <- 26
npop2 <- 30
npop3 <- 34
source("../../R/readms.output.R")
a<- read.ms.output(file="input/msout")
p1.d <- unlist((sapply(a$gam,function(x) colSums(x[1:npop1,]))))
p2.d <- unlist((sapply(a$gam,function(x) colSums(x[(npop1+1):(npop1+npop2),]))))
p3.d <- unlist((sapply(a$gam,function(x) colSums(x[-c(1:(npop1+npop2)),]))))
sfs.p1p2 <- sapply(0:npop1,function(x) table(factor(p2.d[p1.d==x],levels=0:npop2)))
sfs.p1p3 <- sapply(0:npop1,function(x) table(factor(p3.d[p1.d==x],levels=0:npop3)))
sfs.p2p3 <- sapply(0:npop2,function(x) table(factor(p3.d[p2.d==x],levels=0:npop3)))

##2dsfs
write.table(t(sfs.p1p2), 'input/pop1.pop2.true',eol='\t',col.names=F,row.names=F,sep='\t')
write.table(t(sfs.p1p3), 'input/pop1.pop3.true',eol='\t',col.names=F,row.names=F,sep='\t')
write.table(t(sfs.p2p3), 'input/pop2.pop3.true',eol='\t',col.names=F,row.names=F,sep='\t')
##1d
write.table(as.numeric(table(p1.d)), 'input/pop1.true',eol='\t',col.names=F,row.names=F,sep='\t')
write.table(as.numeric(table(p2.d)), 'input/pop2.true',eol='\t',col.names=F,row.names=F,sep='\t')
write.table(as.numeric(table(p3.d)), 'input/pop3.true',eol='\t',col.names=F,row.names=F,sep='\t')


##2dsfs fold
write.table(na2zero(fold2d(t(sfs.p1p2))), 'input/pop1.pop2.fold.true',eol='\t',col.names=F,row.names=F,sep='\t')
write.table(na2zero(fold2d(t(sfs.p1p3))), 'input/pop1.pop3.fold.true',eol='\t',col.names=F,row.names=F,sep='\t')
write.table(na2zero(fold2d(t(sfs.p2p3))), 'input/pop2.pop3.fold.true',eol='\t',col.names=F,row.names=F,sep='\t')


##plot comparisons
if(FALSE){
    pdf("test/fst_folded/log/saf3Tosaf4.folded.fst.comparison.pdf")
    hist((scan('../angsd.master/test/fst_folded/output/pop1.pop2.saf.idx.ml.fold')-scan('test/fst_folded/output/pop1.pop2.saf.idx.ml.fold')),br=101,main="Fold master-dev p1p2")
    hist((scan('../angsd.master/test/fst_folded/output/pop1.pop3.saf.idx.ml.fold')-scan('test/fst_folded/output/pop1.pop3.saf.idx.ml.fold')),br=101,main="Fold master-dev p1p3")
    hist((scan('../angsd.master/test/fst_folded/output/pop2.pop3.saf.idx.ml.fold')-scan('test/fst_folded/output/pop2.pop3.saf.idx.ml.fold')),br=101,main="Fold master-dev p2p3")
    ##
    hist((scan('test/fst_folded/input/pop1.pop2.fold.true')-scan('test/fst_folded/output/pop1.pop2.saf.idx.ml.fold')),br=101,main="Fold true-dev p1p2")
    hist((scan('test/fst_folded/input/pop1.pop2.fold.true')-scan('../angsd.master/test/fst_folded/output/pop1.pop2.saf.idx.ml.fold')),br=101,main="Fold true-master p1p2")
    ##
    hist((scan('test/fst_folded/input/pop1.pop3.fold.true')-scan('test/fst_folded/output/pop1.pop3.saf.idx.ml.fold')),br=101,main="Fold true-dev p1p3")
    hist((scan('test/fst_folded/input/pop1.pop3.fold.true')-scan('../angsd.master/test/fst_folded/output/pop1.pop3.saf.idx.ml.fold')),br=101,main="Fold true-master p1p3")
    ##
    hist((scan('test/fst_folded/input/pop2.pop3.fold.true')-scan('test/fst_folded/output/pop2.pop3.saf.idx.ml.fold')),br=101,main="Fold true-dev p2p3")
    hist((scan('test/fst_folded/input/pop2.pop3.fold.true')-scan('../angsd.master/test/fst_folded/output/pop2.pop3.saf.idx.ml.fold')),br=101,main="Fold true-master p2p3")
    dev.off()
}
