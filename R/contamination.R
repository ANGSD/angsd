controlSNP<-c(-4:-1,1:4)
##this code is horrible. But it works. 
bases=c("A","C","G","T")
require(parallel)

like<-function(x,error,d,freq,eps){
     if(x<0|x>1)
        return(Inf)
    l<-  dbinom(error,d,(1-x)*eps+x*freq)
    return(-sum(log(l)))
}

mom.old <- function(error,d,freq,eps){
    return(mean(error/d-eps)/(mean(freq)-eps))
}

mom.new <- function(error,d,freq,eps){
    return(mean(error/d-eps)/(mean(freq)*(1-4*eps/3)))
}

likeFixed<-function(x,error,d,freq,eps){ 
    if(x<0|x>1)
        return(Inf)
    l<-  dbinom(error,d,x*freq*(1-4*eps/3)+eps)
    return(-sum(log(l)))
}

like1Wrap<-function(x,fixed){
    c<-sum(x[,2])/sum(x[,4])
    r <- c()
    if(fixed==FALSE)
        r <- mom.old(eps=c,error=x[,1],d=x[,3],freq=x[,5])
    else
        r <- mom.new(eps=c,error=x[,1],d=x[,3],freq=x[,5])

    if(fixed==FALSE)
        optimize(like,c(0,1),eps=c,error=x[,1],d=x[,3],freq=x[,5])
    else
        optimize(likeFixed,c(0,1),eps=c,error=x[,1],d=x[,3],freq=x[,5])
}

like1Wrap2<-function(x,max=0.1,fixed){
    c<-sum(x[,2])/sum(x[,4])
    
    if(fixed==FALSE)
        optimize(like,c(0,max),eps=c,error=x[,1],d=x[,3],freq=x[,5])$minimum
    else
        optimize(likeFixed,c(0,max),eps=c,error=x[,1],d=x[,3],freq=x[,5])$minimum
}

jackKnife3<-function(x,fun,mc.cores,...){
   
    ##for matrix per row
    call <- match.call()
    n<-nrow(x) ##n is number of sites
    f<-floor(seq(1,n,n/mc.cores))
    f[mc.cores]<-n+1 ##f contains the start and stop index for the different cores
    ffun<-function(z) unlist(lapply(f[z]:(f[z+1]-1),function(i) fun(x[-i,],...) )) ##we remove index i from each fun call
    u<-unlist(parallel::mclapply(1:(mc.cores-1),ffun,mc.cores=mc.cores))
    thetahat <- fun(x, ...)
    jack.bias <- (n - 1) * (mean(u) - thetahat)
    jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
    obj<-list(jack.se = jack.se, jack.bias = jack.bias, jack.values = u,
              call = call,funhat=thetahat)
    class(obj)<-"blockknife"
    return(obj)
}




estCont<-function(x,jack=FALSE,max=0.1,mc.cores,fixed){
    c<-x$mat[1,2]/sum(x$mat[,2])
  
    err<-sum(x$mat[1,1])/sum(x$mat[,1])
    cat("c est is: ",c," err is",err,"\n")
    ##jx3ack
    j1<-NA
    j3<-NA
    ## jack 2
    
    nS<-length(x$controlSNP)
    
    dat<-cbind(
        err0=matrix(x$error,nrow=nS+1)[rank(c(0,x$controlSNP))[1],],
        err1=colSums(matrix(x$error,nrow=nS+1)[-rank(c(0,x$controlSNP))[1],]),
        d0=matrix(x$d,nrow=nS+1)[rank(c(0,x$controlSNP))[1],],
        d1=colSums(matrix(x$d,nrow=nS+1)[-rank(c(0,x$controlSNP))[1],]),
        freq=x$freq
        )
    
    
    err0=matrix(x$error2,nrow=nS+1)[rank(c(0,x$controlSNP))[1],]
    err1=colSums(matrix(x$error2,nrow=nS+1)[-rank(c(0,x$controlSNP))[1],])
    dat2<-cbind(
        err0,
        err1,
        d0=rep(1,length(err0)),
        d1=rep(nS,length(err1)),
        freq=x$freq
        )
    
    Method1 <- like1Wrap(dat,fixed=fixed)
    Method2 <- like1Wrap(dat2,fixed=fixed)
    
    j12<-NA
    j32<-NA
    if(jack){
        cat ("----------------------\nRunning jackknife for Method1 (could be slow)\n")
        j12<-jackKnife3(dat,like1Wrap2,max=max,mc.cores,fixed=fixed)
        cat ("Running jackknife for Method2 (could be slow)\n")
        j32<-jackKnife3(dat2,like1Wrap2,max=max,mc.cores,fixed=fixed)
    }
    rr <- rbind(cbind(Method1,Method2),se=c(j12[1],j32[1]))
    rownames(rr) <- c("Contamination","llh","SE")
    obj<-list(est=rr,err=err,c=c)
    obj
}

mismatch<-function(r_save,hapMap_save,controlSNP,noNA=TRUE){
    hapPos<-hapMap_save$V1
    cat("mismaatch: nrow(hapPos): ",length(hapPos),"\n")
    pos<-r_save[,1]
   
    keepSites<-pos%in%hapPos
    cat("mismatch nkeeps: ",sum(keepSites),"\n")
    for(c in controlSNP)
        keepSites[keepSites]<-((pos[keepSites]+c)%in%pos)&keepSites[keepSites]
    w<-which(keepSites)
    keep<-sort(c(w,unlist(sapply(w,function(x) x+controlSNP))))
    r<-r_save[keep,]

   
    ##all
    snps<-(r[,1])%in%hapPos
    cat("SNPsites: ",sum(snps),"\n")
    snps1<-!snps
 
    d<-rowSums(r[,-1])#seqdepth
    max<-apply(r[,-1],1,max)#max value
    wmax<-apply(r[,-1],1,which.max) #which base is max
    error<-d-max #max is the difference of non maxobserved
    error2<-rbinom(length(d),1,prob=error/d)
##    Table(error>0)
    keep<-rep(T,nrow(r))
    
    ## test 1
    cat("\n-----------------------\nDoing Fisher exact test for Method1:\n")


    
    mat<-matrix(c(sum(error[snps&keep]),sum(d[snps&keep]-error[snps&keep]),sum(error[snps1&keep]),sum(d[snps1&keep]-error[snps1&keep])),2)

    
    
    colnames(mat)<-c("SNP site","adjacent site")
    rownames(mat)<-c("minor base","major base")
    print(mat)
    print(f<-fisher.test(mat))
    table(wmax[snps])
    table(wmax[snps1])
    
    ## test 2
    cat("\n-----------------------\nDoing Fisher exact test for Method2:\n")
    mat2<-matrix(c(sum(error2[snps&keep]),sum(snps&keep)-sum(error2[snps&keep]),sum(error2[snps1&keep]),sum(snps1&keep)-sum(error2[snps1&keep])),2)

   colnames(mat2)<-c("SNP site","adjacent site")
   rownames(mat2)<-c("minor base","major base")

    print(mat2)
    print(f<-fisher.test(mat2))
    
  
    (mat3<-rbind(
        tapply(error,rep(sort(c(0,controlSNP)),sum(snps)),sum),
        tapply(d-error,rep(sort(c(0,controlSNP)),sum(snps)),sum)
        ))

     cat("\n-----------------------\n major and minor bases - Method1:\n")
    colnames(mat3)[5]<-c("SNP site")
    colnames(mat3)[-5]<-controlSNP

    rownames(mat3)<-c("minor base","major base")
   print(mat3)
     
    (mat4<-rbind(
        tapply(error2,rep(sort(c(0,controlSNP)),sum(snps)),sum),
        tapply(rep(1,length(d))-error2,rep(sort(c(0,controlSNP)),sum(snps)),sum)
        ))
   #  colnames(mat4)<-c("SNP site","adjacent site")[c(2,2,2,2,1,2,2,2,2)]
    colnames(mat4)[5]<-c("SNP site")
    colnames(mat4)[-5]<-controlSNP
    rownames(mat4)<-c("minor base","major base")


     cat("\n-----------------------\n major and minor bases - Method2:\n")
     print(mat4)
    ## get frequencies
    hapMap<-hapMap_save[(hapMap_save[,1]+0 )%in%r[snps,1],] ##<- + Zero because of hg18/hg19
    if(TRUE){
        table(bases[wmax[snps]]==hapMap$V5|bases[wmax[snps]]==hapMap$V2)
        flip<-bases[wmax[snps]]==hapMap$V2
        hapMap$matchAF<-hapMap$V3
        hapMap$matchAF[flip]<-1-hapMap$matchAF[flip]
        freq<-hapMap$matchAF #frequency of the none benny allele
      #  freq[hapMap$V4=="-"]<-1-freq[hapMap$V4=="-"]
    }
##    freq <-1- hapMap[,5]
    obj<-list(error=error,error2=error2,d=d,freq=freq,wmax=wmax,snps=snps,snps1=snps1,mat=mat,mat2=mat2,mat3=mat3,mat4=mat4,controlSNP=controlSNP,r=r)
}


readDat<-function(fileName,maxDepth,minDepth,nSites=1e8){
  to.read<-gzfile(fileName,"rb")
  r_save<-matrix(readBin(to.read,integer(),5*nSites),ncol=5,by=T)
  cat("nSites from ANGSD: ",nrow(r_save),"\n")
  close(to.read)
  dep<-rowSums(r_save[,-1])
  r_save<-r_save[dep>=minDepth&dep<=maxDepth,]
  r_save[,1]<-r_save[,1]-1
  cat("nSites from ANGSD after depfilter: ",nrow(r_save),"\n")
  r_save
}

readHap<-function(MinDist=10,hapFile,minmaf,startPos,stopPos) {
    hapMap_save<-read.table(hapFile,as.is=T)
    cat("HapMap sites:",nrow(hapMap_save), "from file: ",hapFile,"\n")
    hapMap_save <- hapMap_save[!(hapMap_save[,1]<startPos|hapMap_save[,1]>stopPos),]
    cat("HapMap after removing sites outside startPos stopPos ",nrow(hapMap_save),"\n")
    hapMap_save<-hapMap_save[!(hapMap_save[,3]<minmaf|1-hapMap_save[,3]<minmaf),]
    ##    write.table(hapMap_save,file="delme.txt",row.names=F,col.names=F,quote=F)
    cat("HapMap after filtering out minmaf sites",nrow(hapMap_save),"\n")



    hapMap_save <- hapMap_save[(hapMap_save[,2] %in% bases) &(hapMap_save[,5] %in% bases),]
    cat("HapMap after removing undefined(N/-/n) snps ",nrow(hapMap_save),"\n")
    
    hapMap_save <- hapMap_save[!duplicated(hapMap_save[,1]),]
    hapMap_save<-hapMap_save[order(hapMap_save[,1]),]
    cat("unique HapMap sites:",nrow(hapMap_save), "from file: ",hapFile,"\n")
    hapMap_save<-hapMap_save[-which(diff(hapMap_save[,1])<MinDist),]
    ##    write.table(hapMap_save,file="delme.txt",row.names=F,col.names=F,quote=F)
    cat("HapMap after removing close snpts ",nrow(hapMap_save),"\n")

    if(any(is.na(hapMap_save))){
        stop("NA in hapmap")
    }
    #now flip
    fl <- hapMap_save[,4]=="-"
    hapMap_save[hapMap_save[,2]=="A"&fl,2] <- "N"
    hapMap_save[hapMap_save[,2]=="T"&fl,2] <- "A"
    hapMap_save[hapMap_save[,2]=="N",2] <- "T"
    
    hapMap_save[hapMap_save[,2]=="C"&fl,2] <- "N"
    hapMap_save[hapMap_save[,2]=="G"&fl,2] <- "C"
    hapMap_save[hapMap_save[,2]=="N",2] <- "G"
    
    hapMap_save[hapMap_save[,5]=="A"&fl,5] <- "N"
    hapMap_save[hapMap_save[,5]=="T"&fl,5] <- "A"
    hapMap_save[hapMap_save[,5]=="N",5] <- "T"
    
    hapMap_save[hapMap_save[,5]=="C"&fl,5] <- "N"
    hapMap_save[hapMap_save[,5]=="G"&fl,5] <- "C"
    hapMap_save[hapMap_save[,5]=="N",5] <- "G"
    
    return(hapMap_save)
}



if(FALSE){
    mapFile = NULL
    hapFile = NULL
    countFile = NULL
    minDepth=2
    maxDepth=20
    mc.cores=10
    countFile="angsdput.icnts.gz"
    hapFile="RES/HapMapChrX.gz"
    fileName <- "angsdput.icnts.gz"
    fixed=TRUE
    jack=TRUE
    minmaf=0.05
    startPos = 5e6
    stopPos =  154900000
}

doAnal <- function(mapFile,hapFile,countFile,minDepth,maxDepth,mc.cores,fixed,jack,minmaf,startPos,stopPos){
    hapMap_save<-readHap(hapFile=hapFile,minmaf=minmaf,startPos=startPos,stopPos=stopPos)
    r_save<-readDat(countFile,maxDepth,minDepth)

    if(!missing(mapFile)){
        temp<-read.table(mapFile)
        map100<-unlist(apply(temp[,1:2],1,function(x) seq(x[1],x[2],by=1)))
        keep<-r_save[,1]%in%map100
        r_save<-r_save[keep,]
    }
   # maxPos<-154900000
   # minPos<-5e6
    r_save<-r_save[r_save[,1]>=startPos&r_save[,1]<=stopPos,]
    
    res<-mismatch(r_save,hapMap_save,controlSNP)
    res$mat3
    est <-estCont(res,jack=jack,mc.cores=mc.cores,fixed=fixed) 
    print(est$est)
   # est
}



l<-commandArgs(TRUE)
getArgs<-function(x,l)
  unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]

Args<-function(l,args){
    if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
        cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument\n")
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
    cat("->  Needed arguments:\n")
    mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
    cat("->  Optional arguments (defaults):\n")
    mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
    q("no")
}

## choose your parameters and defaults
## NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments

args<-list(
    mapFile = NA,
    hapFile = NULL,
    countFile = NULL,
    minDepth=2,
    maxDepth=20,
    mc.cores=10,
    fixed=TRUE,
    jack=TRUE,
    minmaf=0.05,
    startPos = 5e6,
    stopPos =  154900000,
    seed = NA
        )
##if no argument are given prints the need arguments and the optional ones with default

des<-list(
    mapFile = "Mappability file",
    hapFile = "HapMapfile",
    countFile = "Count file",
    minDepth= "Minimum depth",
    maxDepth= "Maximium depth",
    mc.cores= "Number of cores",
    fixed = "Use fixed version of likelihood",
    jack = "Jacknive to get confidence intervals",
    minmaf = "minimum maf",
    startPos = "start position",
    stopPos = "stop position",
    seed = "set a seed (supply int value)"
       )

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


cat("mapFile = ", mapFile,"\n")
cat("hapFile = ", hapFile,"\n")
cat("countFile =",countFile,"\n")
cat("minDepth = ",minDepth,"\n")
cat("maxDepth = ",maxDepth,"\n")
cat("mc.cores = ",mc.cores,"\n")
cat("fixed = ",fixed,"\n")
cat("jack = ",jack,"\n")
cat("minmaf = ",minmaf,"\n")
cat("startPos = ",startPos,"\n")
cat("stopPos = ",stopPos,"\n")
cat("seed = ",seed,"\n")
{
  if(!is.na(seed))
    set.seed(seed)
    if(!is.na(mapFile))
        doAnal(mapFile=mapFile,hapFile=hapFile,countFile=countFile,minDepth=as.numeric(minDepth),maxDepth=as.numeric(maxDepth),mc.cores=as.numeric(mc.cores),fixed=fixed,jack=jack,minmaf=minmaf,startPos=startPos,stopPos=stopPos)
    else
        doAnal(hapFile=hapFile,countFile=countFile,minDepth=as.numeric(minDepth),maxDepth=as.numeric(maxDepth),mc.cores=as.numeric(mc.cores),fixed=fixed,jack=jack,minmaf=minmaf,startPos=startPos,stopPos=stopPos)
}
