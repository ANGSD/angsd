
MyLike <-function(x,counts){
#    print(x)
    p0011 <- x[1]
    p0111 <- x[2]
    if(p0011+p0111>1)
        return(1e17)
    c1 <- x[3]
    c2 <- x[4]
    
    nAaAa <- counts[1]
    nAaAA <- counts[2]
    nAAaa <- counts[3]
    nAAAa <- counts[4]
    nAAAA <- counts[5]
    en <-   nAaAa*log((2 - 2 *c1 - 2 *c2 + 2 *c1 *c2) *p0011/3)
    to <-   nAaAA*log((4* c2* p0011 - 4 *c1* c2* p0011 + 3 *p0111 - 3 *c1* p0111)/6)
    tre <-  nAAaa*log((4 *p0011 + 4 *c1 *c2 *p0011 + 3 *c1 *p0111 + 3 *c2 *p0111)/12)
    fire<-  nAAAa*log((4* c1* p0011 - 4* c1* c2* p0011 + 3 *p0111 - 3 *c2* p0111)/6)
    fem <-  nAAAA*log((12 - 12 *p0011 + 4* c1* c2* p0011 - 12* p0111 + 3 *c1* p0111 + 3* c2* p0111)/12)
    res <- en+to+tre+fire+fem
    return(-res)
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
    sfsfile = NULL,
    tole=1e-8,
    seed = NA
)
##if no argument are given prints the need arguments and the optional ones with default

des<-list(
    sfsfile = "2sfs from realSFS",
    tole = "tolerance used by optimzation",
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


cat("sfsfile = ", sfsfile,"\n")
cat("tole = ", tole,"\n")
{
    sfs <- read.table(sfsfile)
    if(nrow(sfs)!=1)
        stop("Must supply a single sfs estimate (one line output from realSFS)")
    sfs <- as.numeric(sfs)

    nAaAa <- sfs[5]
    nAaAA <- sfs[4]+sfs[6]
    nAAaa <- sfs[4]+sfs[7]
    nAAAa <- sfs[2]+sfs[8]
    nAAAA <- sfs[1]+sfs[9]
#    nAaAa <- counts[1]
#    nAaAA <- counts[2]
#    nAAaa <- counts[3]
#    nAAAa <- counts[4]
#    nAAAA <- counts[5]

    if(length(sfs)!=9)
        stop("Your SFS must be a joint 2dsfs for two individuals (9 values)")
    cnts <- c(nAaAa,nAaAA,nAAaa,nAAAa,nAAAA)
    print(cnts)
    if(!is.na(seed))
        set.seed(seed)
    ##cnts <- c(37734.380108, 102783.363072, 53170.407313, 160196.133254, 7400237.716253)
    optim(c(.1,.1,.1,.1), MyLike, lower=rep(1e-8,4), upper=rep(1-1e-8,4), method="L-BFGS-B",counts=cnts)
}

#MyLike(counts=c(37734.380108, 102783.363072, 53170.407313, 160196.133254, 7400237.716253),x=c(0.0108979682008016, 0.0379939032592493, 0.309663282543016, 0.0275632180626564))



