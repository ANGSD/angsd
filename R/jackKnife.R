########### do not change ################3
l<-commandArgs(TRUE)
getArgs<-function(x,l)
  unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
Args<-function(l,args){
 if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
  cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
  q("no")
}
 arguments<-list()
 for(a in names(args))
   arguments[[a]]<-getArgs(a,l)

 if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
   cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
   q("no")
 }
 for(a in names(args))
   if(is.null(arguments[[a]]))
     arguments[[a]]<-args[[match(a,names(args))]]

   
 arguments
}

print.args<-function(args,des){
  if(missing(des)){
    des<-as.list(rep("",length(args)))
    names(des)<-names(args)
  }
  cat("->  needed arguments:\n")
  mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("->  optional arguments (defaults):\n")
  mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}
###### ####### ###### ###### ###### #######
# choose your parameters and defaults
# NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(file=NULL,indNames=NULL,outfile="out",boot=0)
#if no argument are given prints the need arguments and the optional ones with default
des<-list(file="the .abbababa filename",indNames="list of individual names (you can use the bam.filelist)",outfile="name of output file",boot="print results for each bootstrap(jackknife), 0=NO")
######################################
#######get arguments and add to workspace
### do not change
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
  cat(" Arguments: output prefix\n")
  q("no")
}
###################################
#file<-"out.abbababa"
#indNames="smallBam.filelist"
#outfile<-"out"
#boot="0"
ab<-read.table(file)
ind<-basename(scan(indNames,what="TheFck"))
N<-length(ind)
boot <- boot!=0

if(N*(N-1)*(N-2)!=ncol(ab)-3){
    cat("The length",N," of indNames does not match the number of individuals in",file,"\n")
    q("n")
}


############functions#################
blockJackUneven<-function(dat){
  nblocks<-nrow(dat)
  X<-cbind(dat[,1]-dat[,2],dat[,1]+dat[,2])
  theta<-function(x) sum(x[,1])/sum(x[,2])
  thetaEst<-theta(X)
  etai<-rep(0, nblocks)
  blockSize<-X[,2]
  blockFrac<-blockSize/sum(blockSize)

  for(i in 1:nblocks)
    etai[i]<-theta(X[-i,])
 
  meanJack<-mean(etai)
  jackEst<-nblocks*thetaEst-sum((1-blockFrac)*etai)
#  jackVar<-(nblocks-1)/nblocks * sum((etai-meanJack)^2)
  jackVar<-1/nblocks * sum( 1/(1/blockFrac-1) * (1/blockFrac*thetaEst-(1/blockFrac-1)*etai - nblocks*thetaEst+sum((1-blockFrac)*etai))^2)
  
 return(c(jackVar=jackVar,jackEst=jackEst,thetaEst=thetaEst))
}
########################################3

Ncomp<-ncol(ab)-3
cat("H1","H2","H3","nABBA","nBABA",c("Dstat","jackEst","SE","Z"), '\n', file=paste(outfile,".txt",sep=""), append=FALSE ,sep="\t")


comp=0
for(h3 in 1:N){
  for(h2 in 1:N){
    if(h2==h3)
      next
    for(h1 in 1:N){
      if(h1==h3)
        next
      if(h1>=h2)
        next
      dat<-ab[,comp*2+3+1:2]
      keep<-rowSums(dat)>0

      if(sum(keep)<3){
          cat("h1 h2 h3 =",h1,h2,h3," has less than 3 blocks. skipping\n")
          jackRes<-c(NA,NA,NA)
          eachD <- c(NA)
      }
      else{

          dat<-dat[keep,]
          jackRes<-blockJackUneven(dat)
          if(boot){
              cc <-colSums(dat)
              eachD <- ((cc[1]-dat[,1]) - (cc[2] - dat[,2]))/ (cc[1]+cc[2] - dat[,1] - dat[,2])
          }              
         
      }
      est<-jackRes[3]
      sd<-sqrt(jackRes[1])
      jack<-jackRes[2]
      Z<-est/sd

      # - print out all estimates
      cat(ind[h1],ind[h2],ind[h3],colSums(dat),c(est,jack,sd,Z), '\n', file=paste(outfile,".txt",sep=""), append=TRUE ,sep="\t")
      if(boot){
          cat(ind[h1],ind[h2],ind[h3],eachD, '\n', file=paste(outfile,".boot",sep=""), append=TRUE ,sep="\t")
      }
      comp<-comp+1
  }}}



cat("output files:\n")
cat("\t",paste(outfile,".txt",sep=""),"\n")
if(boot)
    cat("\t",paste(outfile,".boot",sep=""),"\n")
