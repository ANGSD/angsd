###this code is horrible. But it works. 
bases=c("A","C","G","T")
require(parallel)

like<-function(x,error,d,freq,eps){ 
    if(x<0|x>1)
        return(Inf)
    l<-  dbinom(error,d,(1-x)*eps+x*freq)
    return(-sum(log(l)))
}

like1Wrap2<-function(x,max=0.1){
    c<-sum(x[,2])/sum(x[,4])
    optimize(like,c(0,max),eps=c,error=x[,1],d=x[,3],freq=x[,5])$minimum
}

jackKnife3<-function(x,fun,mc.cores,...){
   
    ##for matrix per row
    call <- match.call()
    n<-nrow(x)
    f<-floor(seq(1,n,n/mc.cores))
    f[mc.cores]<-n+1
    ffun<-function(z) unlist(lapply(f[z]:(f[z+1]-1),function(i) fun(x[-i,],...) ))
    u<-unlist(parallel::mclapply(1:(mc.cores-1),ffun))
    thetahat <- fun(x, ...)
    jack.bias <- (n - 1) * (mean(u) - thetahat)
    jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
    obj<-list(jack.se = jack.se, jack.bias = jack.bias, jack.values = u,
              call = call,funhat=thetahat)
    class(obj)<-"blockknife"
    return(obj)
}



like1Wrap<-function(x){
    c<-sum(x[,2])/sum(x[,4])
    optimize(like,c(0,1),eps=c,error=x[,1],d=x[,3],freq=x[,5])
}

estCont<-function(x,jack=FALSE,max=0.1,mc.cores){
    c<-x$mat[1,2]/sum(x$mat[,2])
    keep<-x$d>0
    err<-sum(x$mat[1,1])/sum(x$mat[,1])

    ##jack
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
    
    Method1 <- like1Wrap(dat)
    Method2 <- like1Wrap(dat2)
    
    j12<-NA
    j32<-NA
    if(jack){
        cat ("----------------------\nRunning jackknife for Method1 (could be slow)\n")
        j12<-jackKnife3(dat,like1Wrap2,max=max,mc.cores)
        cat ("Running jackknife for Method2 (could be slow)\n")
        j32<-jackKnife3(dat2,like1Wrap2,max=max,mc.cores)
    }
    rr <- rbind(cbind(Method1,Method2),se=c(j12[1],j32[1]))
    rownames(rr) <- c("Contamination","llh","SE")
    obj<-list(est=rr,err=err,c=c)
    obj
}

mismatch<-function(r_save,hapMap_save,controlSNP,noNA=TRUE){
    hapPos<-hapMap_save$V2
    
    pos<-r_save[,1]
    {
        if(!noNA){#any
            keepSites<-pos%in%hapPos
            for(c in controlSNP)
                keepSites<-((pos+c)%in%hapPos)|keepSites
            r<-r_save[keepSites,]
        }
        else{
            keepSites<-pos%in%hapPos
            for(c in controlSNP)
                keepSites[keepSites]<-((pos[keepSites]+c)%in%pos)&keepSites[keepSites]
            w<-which(keepSites)
            keep<-sort(c(w,unlist(sapply(w,function(x) x+controlSNP))))
            r<-r_save[keep,]
        }
    }
   
    ##all
    snps<-(r[,1])%in%hapPos
    snps1<-!snps
 
    d<-rowSums(r[,-1])
    max<-pmax(r[,2],r[,3],r[,4],r[,5])
    wmax<-apply(r[,-1],1,which.max)
    error<-d-max
    set.seed(1)
    error2<-rbinom(length(d),1,prob=error/d)
    table(error>0)
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
    hapMap<-hapMap_save[(hapMap_save[,2]+0 )%in%r[snps,1],] ##<- + zero because of hg18/hg19
    table(bases[wmax[snps]]==hapMap$V7|bases[wmax[snps]]==hapMap$V8)
    flip<-bases[wmax[snps]]==hapMap$V8
    hapMap$matchAF<-hapMap$V9
    hapMap$matchAF[flip]<-1-hapMap$matchAF[flip]
    freq<-hapMap$matchAF #frequency of the none benny allele
    freq[hapMap$V6=="-"]<-1-freq[hapMap$V6=="-"]
    obj<-list(error=error,error2=error2,d=d,freq=freq,wmax=wmax,snps=snps,snps1=snps1,mat=mat,mat2=mat2,mat3=mat3,mat4=mat4,controlSNP=controlSNP,r=r)
}


readDat<-function(fileName,maxDepth,minDepth,nSites=1e8){
  to.read<-gzfile(fileName,"rb")
  r_save<-matrix(readBin(to.read,integer(),5*nSites),ncol=5,by=T)
  close(to.read)
  dep<-rowSums(r_save[,-1])
  r_save<-r_save[dep>=minDepth&dep<=maxDepth,]
  r_save[,1]<-r_save[,1]-1
  r_save
}

readHap<-function(MinDist=10,filename) {
    hapMap_save<-read.table(filename)
    hapMap_save<-hapMap_save[order(hapMap_save[,2]),]
    hapMap_save<-hapMap_save[-which(diff(hapMap_save[,2])<MinDist),]
    
    return(hapMap_save)
}

controlSNP<-c(-4:-1,1:4)


##mapFile="../RES/map100.chrX.bz2"
##countFile="/space/anders/ida/idaSjov/kostenkitest/contamination/out/V1countKostinki.USER.bam.X.gz"
##hapFile="../RES/hapMapCeuXlift.map.bz2/"
##fileName <- "angsdput.icnts.gz"

doAnal <- function(mapFile,hapFile,countFile,minDepth,maxDepth,mc.cores){
 
    temp<-read.table(mapFile)
    map100<-unlist(apply(temp[,2:3],1,function(x) seq(x[1],x[2],by=1)))
    r_save<-readDat(countFile,maxDepth,minDepth)

    maxPos<-154900000
    minPos<-5e6

    r_save<-r_save[r_save[,1]>minPos&r_save[,1]<maxPos,]
    keep<-r_save[,1]%in%map100
    r_save<-r_save[keep,]

    
    hapMap_save<-readHap(filename=hapFile)
    res<-mismatch(r_save,hapMap_save,controlSNP)
    res$mat3
    

    est <-estCont(res,jack=T,mc.cores=mc.cores) 
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
    mapFile = NULL,
    hapFile = NULL,
    countFile = NULL,
    minDepth=2,
    maxDepth=20,
    mc.cores=10
    )
##if no argument are given prints the need arguments and the optional ones with default

des<-list(
    mapFile = "Mappability file",
    hapFile = "HapMapfile",
    countFile = "Count file",
    minDepth= "Minimum depth",
    maxDepth= "Maximium depth",
    mc.cores= "Number of cores"
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


doAnal(mapFile=mapFile,hapFile=hapFile,countFile=countFile,minDepth=as.numeric(minDepth),maxDepth=as.numeric(maxDepth),mc.cores=as.numeric(mc.cores))
