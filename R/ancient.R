library("VGAM")
bases <- c("A", "C", "G" ,"T")

collapso<-function(formula,...){
  m<-model.frame(formula)
  mt <- attr(m, "terms")
  x <- model.matrix(mt, m,NULL)
  y<-model.response(m)
  w <- rep(1, nrow(x))

  ty<-.Internal(paste(m[-1], sep="-", collapse=NULL))
  uty<-unique(ty)
  ffun<-function(x){
    if(!is.matrix(x))
      x<-matrix(x,ncol=ncol(y))
    colSums(x)
  }
  yNew<-t(sapply(1:length(uty),function(z) ffun(y[ty==uty[z],])))
  yNew<-matrix(as.integer(yNew),ncol=ncol(y))
  covKeep<-sapply(1:length(uty),function(z) which(ty==uty[z])[1])
  xNew<-x[covKeep,]

  newList<-lapply(m[-1],function(x) x[covKeep])

  attr(x, "assign") = VGAM:::attrassigndefault(x, mt)
  function.name <- "vglm"
  family <- multinomial()
  control = vglm.control(...)
  eval(vcontrol.expression)
  wNew <- rep(1, nrow(xNew))
  attr(xNew, "assign") = VGAM:::attrassigndefault( model.matrix(mt, m,NULL), mt)
  colnames(xNew)<-colnames(x)
  ffit<-vglm.fit(x = xNew, y = yNew, w = wNew,family = multinomial(), control = control,criterion = control$criterion,extra = list(), qr.arg = TRUE, Terms = mt)
  ##ffit<-vglm.fit(x = x, y = y, w = w,family = multinomial(), control = control,criterion = control$criterion,extra = list(), qr.arg = TRUE, Terms = mt)
  ffit$yNew<-yNew
  ffit$newList<-newList
  ffit
}
environment(collapso) <- environment(vglm)

probs<-function(eta){
  ### assumes that est contains one value of eta for each A,C,G
  pT<-1/(1+exp(eta[1])+exp(eta[2])+exp(eta[3]))
  pA<-exp(eta[1])*pT
  pC<-exp(eta[2])*pT
  pG<-exp(eta[3])*pT
  res<-c(pA,pC,pG,pT)
  names(res)<-c("A","C","G","T")
  return(res)
}

funner <- function(){
  pdf("damplots.pdf",width=4*7,height=4*7)
  ps <- min(posi):max(posi)
  par(mfrow=c(4,4))
  bases <- c("A","C","G","T")
  colo <- 1:4
  for(s in c(0,1))
    for(r in 1:4){
      m1 <- t(sapply(ps,function(x) norm(colSums(observed[posi==x&strand==s&ref==r,]))))[,-r]
      m2 <- t(sapply(ps,function(x) norm(colSums(observed[isop==x&strand==s&ref==r,]))))[,-r]
      matplot(ps,m1,type='l',lwd=2,main=paste("ref:",bases[r],"strand:",s,"From beg of read"),col=colo[-r])
      legend("top",bases,fill=colo,box.lty=0)
      matplot(ps,m2,type='l',lwd=2,main=paste("ref:",bases[r],"strand:",s,"From end of read"),col=colo[-r])
      legend("top",bases,fill=colo,box.lty=0)
     # stop("asdfadsafsd")
    }
  dev.off()
}



pplot<-function(x,add=FALSE,turn=FALSE,...){
  reff<-which.max(rowSums(x))
  y<-x[-reff,]
  mmax<-max(y)
  col<-(1:4)[-reff]
  nam<-paste(bases[reff],"->",bases)
  index<-1:ncol(y)
  ylab<-"Error rate"
  if(turn){
    index<-max(index)-index+1
  }
  if(add)
  lines(index,y[1,],col=col[1],...)
    else
  plot(index,y[1,],type="l",ylim=c(0,mmax),col=col[1],ylab=ylab,...)
  lines(index,y[2,],col=col[2],...)
  lines(index,y[3,],col=col[3],...)
legend("top",nam[-reff],lty=1,col=col,bty="n")
}



if(FALSE){

 # source('asdf.R')
  ## Read data and keep only nonzero rows
  thorfinn<-read.table("angsdput.mismatch.gz",header=TRUE)
  thorfinn<-as.matrix(thorfinn)
  
  ref<-thorfinn$Ref
  posi<-thorfinn$posi
  isop<-thorfinn$isop
  qs<-thorfinn$qs
  strand<-thorfinn$strand
  observed <- thorfinn[,-c(1:5)]
  
  ## Only non zero rows
  keepMiss<-rowSums(observed)>0

  ### Reproduce mapdamage plots

  choseRef<-2
  chooseStrand<-0
  keep<-choseRef==ref&chooseStrand==strand&keepMiss
  countsSel<-observed[keep,]
  posiSel<- as.factor(posi[keep])
  isopSel<- as.factor(isop[keep])
  system.time(fit <- collapso( countsSel ~ posiSel-1))

  etas<-coef(fit, matrix=TRUE)
  etas<-matrix(etas, ncol=3,by=TRUE,dimnames=list(sapply(strsplit(names(etas)[1:(length(etas)/3)*3-2],":"),function(x)x[1]),bases[-4]))
### Comparing effects of errors
  len=25
 par(mfrow=c(1,2))
  pplot(sapply(1:len,function(x) probs(etas[x,])),xlab="Posi/isop")

 ##from back

  choseRef<-3
  chooseStrand<-0
  keep<-choseRef==ref&chooseStrand==strand&keepMiss
  countsSel<-observed[keep,]
  posiSel<- as.factor(posi[keep])
  isopSel<- as.factor(isop[keep])
  system.time(fit <- collapso( countsSel ~ isopSel-1))

  etas<-coef(fit, matrix=TRUE)
  etas<-matrix(etas, ncol=3,by=TRUE,dimnames=list(sapply(strsplit(names(etas)[1:(length(etas)/3)*3-2],":"),function(x)x[1]),bases[-4]))
### Comparing effects of errors
  len=25
  pplot(sapply(1:len,function(x) probs(etas[x,])),xlab="Posi/isop")

### Reproduce mapdamage plots

  choseRef<-2
  chooseStrand<-0
  keep<-choseRef==ref&chooseStrand==strand&keepMiss
  countsSel<-observed[keep,]
  posiSel<- as.factor(posi[keep])
  isopSel<- as.factor(isop[keep])
  system.time(fit <- collapso( countsSel ~ posiSel-1))

  etas<-coef(fit, matrix=TRUE)
  etas<-matrix(etas, ncol=3,by=TRUE,dimnames=list(sapply(strsplit(names(etas)[1:(length(etas)/3)*3-2],":"),function(x)x[1]),bases[-4]))
### Comparing effects of errors
  len=25
 par(mfrow=c(1,2))
  pplot(sapply(1:len,function(x) probs(etas[x,])),xlab="Posi/isop")

 
 choseRef<-3
  chooseStrand<-1
  keep<-choseRef==ref&chooseStrand==strand&keepMiss
  countsSel<-observed[keep,]
  posiSel<- as.factor(posi[keep])
  isopSel<- as.factor(isop[keep])
  system.time(fit <- collapso( countsSel ~ isopSel-1))

  etas<-coef(fit, matrix=TRUE)
  etas<-matrix(etas, ncol=3,by=TRUE,dimnames=list(sapply(strsplit(names(etas)[1:(length(etas)/3)*3-2],":"),function(x)x[1]),bases[-4]))
### Comparing effects of errors
  len=25
  pplot(sapply(1:len,function(x) probs(etas[x,])),xlab="Posi/isop",lty=2,add=T)

 #funky stuff below

 
  choseRef<-3
  chooseStrand<-1
  keep<-choseRef==ref&chooseStrand==strand&keepMiss
  countsSel<-observed[keep,]
  posiSel<- as.factor(posi[keep])
  isopSel<- as.factor(isop[keep])
  system.time(fit <- collapso( countsSel ~ posiSel+isopSel))

  etas<-coef(fit, matrix=TRUE)
  etas<-matrix(etas, ncol=3,by=TRUE,dimnames=list(sapply(strsplit(names(etas)[1:(length(etas)/3)*3-2],":"),function(x)x[1]),bases[-4]))
### Comparing effects of errors
 
 ii<-length(unique(isop))
  pp<-length(unique(posi))
  len=21

  pplot(sapply(1:len,function(x) probs(etas[1,]+etas[x+1,]+etas[pp+len-x+1,])),xlab="Posi/isop")

  choseRef<-2
  chooseStrand<-1
  keep<-choseRef==ref&chooseStrand==strand&keepMiss
  countsSel<-observed[keep,]
  posiSel<- as.factor(posi[keep])
  isopSel<- as.factor(isop[keep])
  system.time(fit2 <- collapso( countsSel ~ posiSel+isopSel))

  etas<-coef(fit2, matrix=TRUE)
  etas<-matrix(etas, ncol=3,by=TRUE,dimnames=list(sapply(strsplit(names(etas)[1:(length(etas)/3)*3-2],":"),function(x)x[1]),bases[-4]))
### Comparing effects of errors
  len=25
 ii<-length(unique(isop))
  pp<-length(unique(posi))
  len=50

  pplot(sapply(1:len,function(x) probs(etas[1,]+etas[x+1,]+etas[pp+len-x+1,])),xlab="Posi/isop",add=T)

}
