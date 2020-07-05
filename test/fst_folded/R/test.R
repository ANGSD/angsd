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
