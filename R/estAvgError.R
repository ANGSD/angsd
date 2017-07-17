options(warn=-1)
bases<-c("A","C","G","T")
b <- c(bases,"N")
library(pracma)
library(data.table)

########### do not change ################
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
  cat("->  Needed arguments:\n")
  mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("->  Optional arguments (defaults):\n")
  mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}

## choose your parameters and defaults
## NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(angsdFile = "out",
           errFile = FALSE,
           nameFile = FALSE,
           sizeFile = FALSE,
           out="out",
           addErr=FALSE,
           nIter=100,
           main="",
           maxErr=0.02
           )

#if no argument aree given prints the need arguments and the optional ones with default
des<-list(angsdFile="output angsdFile (no .abbababa2 extension) from 'angsd -doAbbababa2 1' command",
          errFile="list of files of populations (with extension. For groups with no error file write NA in the related line.)",
          nameFile="file with population names",
          sizeFile="file with population sizes",
          out="Name of the out files",
          addErr="amount of error correction to add with increment and bases for transitions. E.g. -0.005,0.005,0.001;A,C;G,T",
          nIter="Number of optimazation attemps",
          maxErr="maximum allowed error rate",
          main="Title for the plots"
          )

####### get arguments and add to workspace
###     do not change
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
  cat(" Arguments: output prefix\n")
  q("no")
}

###################################

cat("----------\nangsdFile: ",angsdFile,"  errFile: ", errFile," nameFile: ",nameFile," sizeFile: ",sizeFile," out: ", out," addErr: ", addErr," nIter: ",nIter,"maxErr: ",maxErr,"\n-----------\n" )

b<-c("A","C","G","T","N")

nIter<-as.integer(nIter)
maxErr=as.numeric(maxErr)

readTable <- function(file){   
    r <- as.matrix(read.table(file, header=FALSE))
    return(r)
}
    
getMat<-function(x){
    m<-array(0,dim=c(5,5,5),dimnames=list(b,b,b))
    for(s in 0:4)
        for(p in 0:4)
            for(a in 0:4)
                m[a+1,p+1,s+1]<-x[a*25+p*5+s+1]
    return(m)
}

logLike<-function(x,Xch,Pch){
    eMat<-matrix(0,4,4)
    eMat[-c(1,6,11,16)]<-x
    diag(eMat)<-1-rowSums(eMat)
    P <- Pch %*% eMat
    ll <- -sum(log(P)*Xch)
    return(ll)
}

getFromErrFile <- function(r,res,maxErr,nInd,logLike){
                                        #options(warn=-1)#suppress warnings
    #print(dim(r))
    r = t(r)
    for(j in 1:nInd){
        m<-getMat(r[j,])
        
        ##remove if missing
        m<-m[-5,-5,-5]
        
        Pch<-matrix(0,4,4)
        for(i in 1:4)
            Pch<-Pch+m[,,i]
        
        Pch<-Pch/rowSums(Pch)
        
        Xch<-matrix(0,4,4)
        for(i in 1:4)
            Xch<-Xch+m[,i,]
        
        conv <- nlminb(runif(12)/100,logLike,upper=rep(maxErr,12),lower=rep(1e-10,12),Xch=Xch,Pch=Pch)
        for(i in 1:nIter){
            Tempconv <- nlminb(runif(12)/100,logLike,upper=rep(maxErr,12),lower=rep(1e-10,12),Xch=Xch,Pch=Pch)
            if(Tempconv$objective<conv$objective){
                conv<-Tempconv               
            }
        }
        #options(warn=0)#turn on warnings again
        return(rbind(res,conv$par))
    }
}

buildMat <- function(res){
    resMat <- matrix(nrow=4,ncol=4)
    resMat[2:4,1] <- res[1:3]
    resMat[c(1,3,4),2] <- res[4:6]
    resMat[c(1,2,4),3] <- res[7:9]
    resMat[1:3,4] <- res[10:12]
    diag(resMat) <- c( 1-sum(res[1:3]),1-sum(res[4:6]),1-sum(res[7:9]),1-sum(res[10:12]) )
    return(resMat)
}



buildInv <- function(res){
    return(solve(res))
}


getErrMat <- function(resMat){
    errMat = matrix(nrow=256,ncol=256)

    idRow = list()
    idCol = list()
    cont = 1
    
    for(i in 1:4){
        for(j in 1:4){
            for(k in 1:4){
                for(l in 1:4){
                    idRow[[cont]] = c(i,j,k,l)
                    idCol[[cont]] = c(i,j,k,l)
                    cont = cont +1}}}}
    
    for(rw in 1:256)
        for(cl in 1:256)
            errMat[rw,cl] = resMat[[1]][idRow[[rw]][1],idCol[[cl]][1]] * resMat[[2]][idRow[[rw]][2],idCol[[cl]][2]] * resMat[[3]][idRow[[rw]][3],idCol[[cl]][3]] * resMat[[4]][idRow[[rw]][4],idCol[[cl]][4]]
    
    return(errMat)
}

#ABBA<-c("X0110","X0220","X0330","X1001","X1221","X1331","X2002","X2112","X2332","X3003","X3113","X3223")
#BABA<-c("X0101","X0202","X0303","X1010","X1212","X1313","X2020","X2121","X2323","X3030","X3131","X3232")
#ABBAtr<-c("X0110","X0330","X1001","X1221","X2112","X2332","X3003","X3223")
#BABAtr<-c("X0101","X0303","X1010","X1212","X2121","X2323","X3030","X3232")
#BBAA<-c("X0011","X0022","X0033","X1100","X1122","X1133","X2200","X2211","X2233","X3300","X3311","X3322")


ABBA<-c("0110","0220","0330","1001","1221","1331","2002","2112","2332","3003","3113","3223")
BABA<-c("0101","0202","0303","1010","1212","1313","2020","2121","2323","3030","3131","3232")
ABBAtr<-c("0110","0330","1001","1221","2112","2332","3003","3223")
BABAtr<-c("0101","0303","1010","1212","2121","2323","3030","3232")
BBAA<-c("0011","0022","0033","1100","1122","1133","2200","2211","2233","3300","3311","3322")


getJackKnife <- function(outData,finalInv=FALSE,ABBAname,BABAname,BBAAname){
    if(sum(as.numeric(outData[,5])==0) == dim(outData)[1])
        return(list(thetaN=NA,thetaJack=NA,varJack=NA,Z=NA,pv=NA,nABBA=NA,nBABA=NA,nBBAA=NA,numBlock=NA))  
    colWeights = as.vector(as.numeric(rowSums(outData[,-c(1:6)])))
    #zeroIdx <- outData[,6]!=0 
    zeroIdx <- colWeights!=0
    outData <- outData[zeroIdx,]
    colWeights <- colWeights[zeroIdx]
    
    skipData=FALSE
    
    if(is.vector(outData)){
        skipData = TRUE} else if(nrow(outData)==0){
        skipData = TRUE}
    
    if(skipData)
	return(list(thetaN=NA,thetaJack=NA,varJack=NA,Z=NA,pv=NA,nABBA=NA,nBABA=NA,nBBAA=NA,numBlock=NA))

    #zeroIdx <- outData[,5]!=0
    colDen <- rowSums(outData[,c(ABBAname)]) + rowSums(outData[,c(BABAname)])
    zeroIdx <- colDen!=0
    outData <- outData[zeroIdx,]
    colWeights <- colWeights[zeroIdx]
    
    if(is.vector(outData)){
        skipData = TRUE}  else if(nrow(outData)==0){
        skipData = TRUE}

    if(skipData)
	return(list(thetaN=NA,thetaJack=NA,varJack=NA,Z=NA,pv=NA,nABBA=NA,nBABA=NA,nBBAA=NA,numBlock=NA))
 
    #seenSites = sum(as.numeric(outData[,6]))
    #weigth = colWeights
    weigth = as.numeric( rowSums(outData[,c(ABBAname,BABAname)])  )
    outData = outData[,-c(1,2,3,4,5,6)]
    Edata <- numeric(prod(dim(outData)))
    dim(Edata)<- dim(outData)
    L = nrow(outData)
    totAbba = 0
    totBaba = 0
    totBbaa = 0
    num = rep(0,L)
    den = rep(0,L)

    if(!is.matrix(finalInv)){
        num = rowSums(outData[,ABBAname]) - rowSums(outData[,BABAname])
        den = rowSums(outData[,ABBAname]) + rowSums(outData[,BABAname])
        totAbba = sum(outData[,ABBAname])
        totBaba = sum(outData[,BABAname])
        totBbaa = sum(outData[,BBAAname])
    } else{
        res = apply(outData,1,
                    function(x){
                        x = as.numeric(x)
                        Edata = finalInv %*% x
                        x = Edata
                        names(x) = colnames(outData)
                        num = sum(as.numeric(x[ABBAname])) - sum(as.numeric(x[BABAname]))
                        den = sum(as.numeric(x[ABBAname])) + sum(as.numeric(x[BABAname]))
                        totAbba =  sum(as.numeric(x[ABBAname]))
                        totBaba =  sum(as.numeric(x[BABAname]))
                        totBbaa =  sum(as.numeric(x[BBAAname]))
                        return(c(num,den,totAbba,totBaba,totBbaa))
                    })
        num = res[1,]
        den = res[2,]
        totAbba = sum(res[3,])
        totBaba = sum(res[4,])
        totBbaa = sum(res[5,])
    }
    
    #block jack knife estimation
    remIdx = (is.nan(num) | is.nan(den)) | (is.na(num) | is.na(den))
    
    num = num[!remIdx]
    den = den[!remIdx]
    weigth = weigth[!remIdx]/sum(weigth[!remIdx]) #block weights
    L = length(num)
    thetaN = sum(num)/sum(den) #D statistics calculated without jackknife
    errorCount = 0
     
    thetaJStar <- rep(0,L) #partial estimators
    thetaJack <- rep(0,L)  #Block Jack knife estimator
    pseudo <- 1/weigth
    thetaJStar <- sapply( 1:L, function(x) sum(as.numeric(num[-x]))/sum(as.numeric(den[-x])) )  #D statistic in each block
    thetaJack <- L*thetaN - sum((1-weigth) * thetaJStar) #jackknife D statistic       
    meanJack <- mean(thetaJStar)# average D statistics over the blocks
    thetaTilde <- pseudo*thetaN-(pseudo-1) * thetaJStar #intermediate quantity for the variance of the jackknife D statistic
    varJack <- 1/L * sum( 1/(pseudo-1) * (thetaTilde - thetaJack)^2 ) #variance of the jackknife D   
    Z <- thetaN / sqrt(varJack) #Z value for standard normal
    pv = 2*min(pnorm(Z,0,1),1-pnorm(Z,0,1)) #pvalue for standard normal
    return(list(thetaN=thetaN,thetaJack=thetaJack,varJack=varJack,Z=Z,pv=pv,nABBA=totAbba,nBABA=totBaba,nBBAA=totBbaa,numBlock=L))
}

angsdFile = paste(angsdFile,".abbababa2",sep="")
#outDataTotal <- read.table(angsdFile,header=T,as.is=T)#old read data
#outDataTotal <- as.matrix(fread(input=angsdFile,sep="\t",showProgress=TRUE,header=TRUE,data.table=TRUE))
erCor= (errFile != FALSE)

if(sizeFile==FALSE && nameFile==FALSE){
	cat("Error: Define at least one between nameFile and sizeFile!\n")
	q("no")
}	


if(sizeFile==FALSE){
	namePop = unlist(read.table(nameFile, header=FALSE))
	L = length(namePop)
	sizePop = rep(1,L)
} else{	
	sizePop = unlist(read.table(sizeFile, header=FALSE, as.is=TRUE))
}

if(nameFile==FALSE){
	sizePop = unlist(read.table(sizeFile, header=FALSE, as.is=TRUE))
	L = length(sizePop)
	namePop = c()
	for(i in 1:L)
		namePop = c(namePop,paste("Population",i,sep="_"))
} else{	
	namePop = unlist(read.table(nameFile, header=FALSE))
}

numPop = length(namePop)
numComb = choose(numPop-2,2)*(numPop-1)
cumPop = c(0,cumsum(sizePop))

#lenList = length(outDataTotal[,1])

addErrors=(addErr!=FALSE)

if(addErrors){
    addErrStr = unlist(strsplit(addErr, split=";"))
    lenStr = length(addErrStr)
    addErrVal = eval(parse(text=paste("c(",addErrStr[1],")",sep="")))
    if(strcmp("",addErrStr[2])){
        addFrom=c(1,2,3,4)
        } else{
            lettersFrom=unlist(strsplit(addErrStr[2], split=","))
            addFrom=c()
            letters = c("A","C","G","T")
            for(l in lettersFrom)
                for(i in 1:4)
                    if(strcmp(l,letters[i]))
                        addFrom=c(addFrom,i)
        }
                
    if(strcmp("",addErrStr[2])){
        addTo=c(1,2,3,4)
        } else{
            lettersTo=unlist(strsplit(addErrStr[3], split=","))
            addTo=c()
            letters = c("A","C","G","T")
            for(l in lettersTo)
                for(i in 1:4)
                    if(strcmp(l,letters[i]))
                        addTo=c(addTo,i)
        }
}

combs = matrix(nrow=numComb, ncol=4)
sizes = matrix(0, nrow=numComb, ncol=4)
cumId = matrix(nrow=numComb, ncol=4)
nameId = matrix(nrow=numComb, ncol=4)
combNames = c()
posit = 1
for(i in 1:(numPop-1)){
    for(j in 1:(numPop-1)){
        for(k in 1:(numPop-1)){
                if(j>i && j!=k && i!=k){
                    combs[posit,1]=i; combs[posit,2]=j; combs[posit,3]=k; combs[posit,4]=numPop;
                    sizes[posit,1]=sizePop[i]; sizes[posit,2]=sizePop[j]; sizes[posit,3]=sizePop[k]; sizes[posit,4]=sizePop[numPop];
                    cumId[posit,1]=cumPop[i]; cumId[posit,2]=cumPop[j]; cumId[posit,3]=cumPop[k]; cumId[posit,4]=cumPop[numPop];
                    combNames=c(combNames, paste("",namePop[i],namePop[j],namePop[k],namePop[numPop],sep="."))
                    nameId[posit,] = c(as.character(namePop[i]),as.character(namePop[j]),as.character(namePop[k]),as.character(namePop[numPop]))
                    posit = posit + 1
            }
        }
    }
}


#get error matrices
if(erCor==1){
    errFiles = unlist(read.table(errFile, header=FALSE, as.is=TRUE))
    cont=1
    resMat = list();
    cat("Estimating error matrix from file:\n")
    for(str in errFiles){    
        nInd<-1
        res=NULL
        resMat[[cont]] = diag(c(1,1,1,1));
        filez=FALSE;     
        if(!is.na(str)){
            resMat[[cont]] = diag(c(0,0,0,0));
            r = readTable(str)
            for( n in 1:nrow(r) ){
                res=NULL
                res = getFromErrFile(as.matrix(r[n,]),res,maxErr,nInd,logLike)
                resMat[[cont]] = resMat[[cont]] + buildMat(res)
                #print(resMat[[cont]])
            }
            resMat[[cont]] = resMat[[cont]] / nrow(r)
            cat(sprintf("\t\"%s\"\n",str))
        }
        cont = cont + 1
    }
}

if(erCor==1){
    for(ii in 1:numPop){
        fileOut = paste(out,".barPlotErrors.",namePop[ ii ],".pdf",sep="")
        pdf(fileOut)
        errors=as.vector(resMat[[ ii ]])[-c(1,6,11,16)]
        
        tickz=max(errors)
        xLabel=c("A-->C","A-->G","A-->T","C-->A","C-->G","C-->T","G-->A","G-->C","G-->T","T-->A","T-->C","T-->G")
        
        par(mar=c(5,5,4.1,2.1))
        barplot(errors,beside=T,col=c("darkgreen"),border=NA,width=9,xlim=c(0,140),ylim=c(0,1.1*tickz),yaxt="n",main=sprintf("Type-specific errors: %s",namePop[ ii ]))
        
        mtext(side=2,line=4,"Error")
        mtext(side=1,line=4,"Type-specific Error")
        axis(side=1, at=seq(10,140,11)-6, labels=xLabel, las=2, cex=0.5, tck=0)

        yTickz=seq(0,1.1*tickz,1.1*tickz/20)
        len=length(yTickz)
        yLabel=c()
        for(i in 1:len)
            yLabel[i]=sprintf("%.1e",yTickz[i])
        yLabel=toString(yLabel)
        yLabel = eval(parse(text=paste("c(",yLabel,")",sep="")))
        yLabel = strtrim(yLabel,7)
        abline(h=seq(0,1.1*tickz,1.1*tickz/20), lwd=0.5, col="gray" )
        axis(side=2, at=seq(0,1.1*tickz,1.1*tickz/20), labels=yLabel, las=2, cex=0.3, tck=0)
        
        dev.off()
    }
}


cat("--- Building tree error matrices ---\n ")
ptm <- proc.time() #start timer
if(erCor==1){
    solveMat = list()
    bigMat = list()
    pb <- txtProgressBar(min=1,max=numComb,initial = 0,style = 3,char=":)", width=20)
    bigMat = list()
    bigMat[[4]] = resMat[[ numPop ]]
    for(idComb in 1:numComb){
        id = combs[idComb,]
        bigMat[[1]] = resMat[[ id[1] ]]
        bigMat[[2]] = resMat[[ id[2] ]]
        bigMat[[3]] = resMat[[ id[3] ]]
        solveMat[[idComb]] =  buildInv( getErrMat(bigMat) )
        setTxtProgressBar(pb,idComb,label="ciaooo")
        flush.console()
    }
}
cat("--- Time Spent ",(proc.time() - ptm)[3]," sec \n")
finalptm = (proc.time() - ptm)[3]

ptm2 <- proc.time() #start timer
FILEOBS<-paste(out,".Observed",".txt",sep="")
FILEERROR<-paste(out,".ErrorCorr",".txt",sep="")
FILETRANS<-paste(out,".TransRem",".txt",sep="")
FILEERRORTRANS<-paste(out,".ErrorCorr.TransRem",".txt",sep="")

outDataTotal <- as.matrix(fread(input=angsdFile,sep="\t",showProgress=TRUE,header=TRUE,data.table=TRUE))
lenList = length(outDataTotal[,1])

for(idComb in 1:numComb){
    
    idxList = seq(1,lenList,numComb) + (idComb-1)
    outData = outDataTotal[idxList,]
    dims <- dim(outData)
    outData <- as.numeric(outData)
    dim(outData) <- dims
    colnames(outData) <- colnames(outDataTotal)
    
    
    id = combs[idComb,]
    sz = sizes[idComb,]
    nm = nameId[idComb,]
            
    if(erCor==1){#ERROR CORRECTED D
        
        result1 = getJackKnife(outData,solveMat[[idComb]],ABBAname=ABBA,BABAname=BABA,BBAAname=BBAA)

        if(idComb==1){
            str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBlocks\tH1\tH2\tH3\tH4")
            str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result1$thetaN,result1$thetaJack,result1$varJack,result1$Z,result1$pv,result1$nABBA,result1$nBABA,result1$numBlock,nm[1],nm[2],nm[3],nm[4])
            cat(str,str2,file=FILEERROR,sep="\n")} else{
                                                     str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result1$thetaN,result1$thetaJack,result1$varJack,result1$Z,result1$pv,result1$nABBA,result1$nBABA,result1$numBlock,nm[1],nm[2],nm[3],nm[4])
                                                     cat(str2,file=FILEERROR,sep="\n",append=TRUE)
            }

        result3 = getJackKnife(outData,solveMat[[idComb]],ABBAname=ABBAtr,BABAname=BABAtr,BBAAname=BBAA)

        if(idComb==1){
            str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBlocks\tH1\tH2\tH3\tH4")
            str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result3$thetaN,result3$thetaJack,result3$varJack,result3$Z,result3$pv,result3$nABBA,result3$nBABA,result3$numBlock,nm[1],nm[2],nm[3],nm[4])
            cat(str,str2,file=FILEERRORTRANS,sep="\n")} else{
                str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result3$thetaN,result3$thetaJack,result3$varJack,result3$Z,result3$pv,result3$nABBA,result3$nBABA,result3$numBlock,nm[1],nm[2],nm[3],nm[4])
                cat(str2,file=FILEERRORTRANS,sep="\n",append=TRUE)
            }
    }

### write in files the effect of added error to transition FROM --> TO
    if(addErrors){###########
        
        addErr = seq(addErrVal[1], addErrVal[2], addErrVal[3])
        LErr = length(addErr)
        letters = c("A","C","G","T")
        Ltot = LErr*length(addTo)*length(addFrom)*3
        idBar=1
        
        if(idComb==1){
            dirName = paste(out,".errorDataFolder",sep="")
            dir.create(dirName, showWarnings = FALSE, recursive = FALSE, mode = "0777")    
        }
        cat("Adding errors to H1, H2 or H3\n")
        pb <- txtProgressBar(min=1,max=Ltot,initial = 0,style = 3,char=":)", width=20)
        for(FROM in addFrom){
            for(TO in addTo){
                if(FROM!=TO){
                    for(ii in 1:3){
                        newMat = list()
                        newMat[[1]]=resMat[[id[1]]]
                        newMat[[2]]=resMat[[id[2]]]
                        newMat[[3]]=resMat[[id[3]]]
                        newMat[[4]]=resMat[[4]]
                        for(er in 1:LErr){
                            newMat[[ii]]=resMat[[id[ii]]]
                            fileOut1 = paste("Add",addErr[er],".H",ii,".",letters[FROM],"2",letters[TO],combNames[ idComb ],".ErrCorr",".txt",sep="")
                            fileOut1 = paste(dirName,"/",fileOut1,sep="")
                            newMat[[ii]][TO,FROM] = newMat[[ii]][TO,FROM] + addErr[er]
                            newMat[[ii]][FROM,FROM] = newMat[[ii]][FROM,FROM] - addErr[er]
                            newErrMat = getErrMat(newMat)
                            newInv = buildInv(newErrMat)
                            result = getJackKnife(outData,newInv,ABBAname=ABBA,BABAname=BABA,BBAAname=BBAA)
                            
                            str = sprintf("NumInd\taddErr\tD\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tH1\tH2\tH3\tH4")
                            str2 = sprintf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%s",ii,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv,result$nABBA,result$nBABA,nm[1],nm[2],nm[3],nm[4])
                            cat(str,str2,file=fileOut1,sep="\n")
                            setTxtProgressBar(pb, idBar,label="ciao"); idBar=idBar+1;
                            gc()
                            flush.console()
                        }
                    }
                }  
            }
        }
        close(pb); idBar=0;
        
        cat("Add errors to H1,H2 or H3 + Error Corr + Trans Removal + Plots\n")
        pb <- txtProgressBar(min=1,max=Ltot,initial = 0,style = 3,char=":)", width=20)
        for(FROM in addFrom){
            for(TO in addTo){
                if(FROM!=TO){
                    for(ii in 1:3){
                        newMat = list()
                        newMat[[1]]=resMat[[id[1]]]
                        newMat[[2]]=resMat[[id[2]]]
                        newMat[[3]]=resMat[[id[3]]]
                        newMat[[4]]=resMat[[4]]
                        for(er in 1:LErr){
                            newMat[[ii]]=resMat[[id[ii]]]
                            fileOut2 = paste("Add",addErr[er],".H",ii,".",letters[FROM],"2",letters[TO],combNames[ idComb ],".ErrCorr.TransRem",".txt",sep="")
                            fileOut2 = paste(dirName,"/",fileOut2,sep="")
                            newMat[[ii]][TO,FROM] = newMat[[ii]][TO,FROM] + addErr[er]
                            newMat[[ii]][FROM,FROM] = newMat[[ii]][FROM,FROM] - addErr[er]
                            newErrMat = getErrMat(newMat)
                            newInv = buildInv(newErrMat)
                            result = getJackKnife(outData,newInv,ABBAname=ABBAtr,BABAname=BABAtr,BBAAname=BBAA)
                            str = sprintf("NumInd\taddErr\tD\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tH1\tH2\tH3\tH4")
                            str2 = sprintf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%s",id[ii],addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv,result$nABBA,result$nBABA,nm[1],nm[2],nm[3],nm[4])
                            cat(str,str2,file=fileOut2,sep="\n")
                            setTxtProgressBar(pb, idBar,label="ciao"); idBar=idBar+1;
                            gc()
                            flush.console()
                        }
                    }
                }
            }     
        }
        close(pb); idBar=0;

        # Make the plots
        # Read the data from output files and plot it
        dat=matrix(NA,nrow=3,ncol=LErr)
        datNT=matrix(NA,nrow=3,ncol=LErr)
        for(FROM in addFrom){
            for(TO in addTo){
                if(FROM!=TO){
                for(ii in 1:3){
                    for(er in 1:LErr){                                      
                        #read files
                        pFile = paste("Add",addErr[er],".H",ii,".",letters[FROM],"2",letters[TO],combNames[ idComb ],".ErrCorr",".txt",sep="")
                        pFile = paste(dirName,"/",pFile,sep="")
                        pFileNT = paste("Add",addErr[er],".H",ii,".",letters[FROM],"2",letters[TO],combNames[ idComb ],".ErrCorr.TransRem",".txt",sep="")
                        pFileNT = paste(dirName,"/",pFileNT,sep="")
                        dat[ii,er] = as.numeric(read.table(pFile,as.is=T,header=T)[3])
                        datNT[ii,er] = as.numeric(read.table(pFileNT,as.is=T,header=T)[3])
                    }
                }
                                        #plot
                yLimits=c( c(min(as.vector(dat),as.vector(datNT))),  c(max(as.vector(dat),as.vector(datNT))))
                xAxis=seq(addErrVal[1],addErrVal[2],addErrVal[3])
                str=sprintf("Effect on D removing transition error on %s-->%s",letters[FROM],letters[TO])

                fileOut = paste("plotAddErr",".",letters[FROM],"2",letters[TO],combNames[ idComb ],".pdf",sep="")
                fileOut = paste(dirName,"/",fileOut,sep="")

                pdf(fileOut)
                plot(c(addErrVal[1],addErrVal[2]),c(0,0),type="l",col = "gray",ylim=yLimits, xlab="Error removed",ylab="D statistic",main=str,lwd=0.5)
                lines(xAxis, dat[1,], type="l",col="darkgreen",lwd=2.5)
                lines(xAxis, dat[2,], type="l",col="darkorange1",lwd=2.5)
                lines(xAxis, dat[3,], type="l",col="darkmagenta",lwd=2.5)
                lines(xAxis, datNT[1,], type="l",col="darkgreen",lwd=2.5,lty=2)
                lines(xAxis, datNT[2,], type="l",col="darkorange1",lwd=2.5,lty=2)
                lines(xAxis, datNT[3,], type="l",col="darkmagenta",lwd=2.5,lty=2)
                legend("top",c(nm[1],nm[2],nm[3],"Err.Corr.","Err.Corr. + Trans.Rem."), lty = c(1,1,1,1,2), lwd=c(2,2,2,2,2), col = c("darkgreen","darkorange1","darkmagenta","black","black"),bty="n",ncol=2)
                dev.off()
                }
            }
        }        



#REMOVE FILES
        for(FROM in addFrom){
            for(TO in addTo){
                if(FROM!=TO){
                    for(ii in 1:3){
                        for(er in 1:LErr){
                            fileOut1 = paste("Add",addErr[er],".H",ii,".",letters[FROM],"2",letters[TO],combNames[ idComb ],".ErrCorr",".txt",sep="")
                            fileOut1 = paste(dirName,"/",fileOut,sep="")
                            file.remove(fileOut1)
                        }
                    }
                }  
            }
        }
        close(pb); idBar=0;
        
        for(FROM in addFrom){
            for(TO in addTo){
                if(FROM!=TO){
                    for(ii in 1:3){
                        for(er in 1:LErr){
                            fileOut2 = paste("Add",addErr[er],".H",ii,".",letters[FROM],"2",letters[TO],combNames[ idComb ],".ErrCorr.TransRem",".txt",sep="")
                            fileOut2 = paste(dirName,"/",fileOut,sep="")
                            file.remove(fileOut2)
                        }
                    }
                }
            }     
        }





        

###bar plots of errors
        if(FALSE){
        if(erCor==1){
            fileOut = paste(dirName,"/barPlotErrors",combNames[ idComb ],".pdf",sep="")
            pdf(fileOut)
            err1=as.vector(bigMat[[1]])[-c(1,6,11,16)]
            err2=as.vector(bigMat[[2]])[-c(1,6,11,16)]
            err3=as.vector(bigMat[[3]])[-c(1,6,11,16)]
            err4=as.vector(bigMat[[4]])[-c(1,6,11,16)]

            tickz=max(err1,err2,err3,err4)
            xLabel=c("A-->C","A-->G","A-->T","C-->A","C-->G","C-->T","G-->A","G-->C","G-->T","T-->A","T-->C","T-->G")
            errors = rbind(err1,err2,err3,err4)

            par(mar=c(5,5,4.1,2.1))
            barplot(errors,beside=T,col=c("darkgreen","darkorange1","darkmagenta","firebrick1"),border=NA,width=2.55,xlim=c(0,160),ylim=c(0,1.2*tickz),yaxt="n",main="Type-specific error rates for the individuals in H1, H2, H3, H4")
            mtext(side=2,line=4,"Error")
            mtext(side=1,line=4,"Type-specific Error")
            axis(side=1, at=seq(10,160,13)-6, labels=xLabel, las=2, cex=0.5, tck=0)

            yTickz=seq(0,1.2*tickz,1.2*tickz/20)
            len=length(yTickz)
            yLabel=c()
            for(i in 1:len)
                yLabel[i]=sprintf("%.1e",yTickz[i])
            yLabel=toString(yLabel)
            yLabel = eval(parse(text=paste("c(",yLabel,")",sep="")))
            yLabel = strtrim(yLabel,7)
            abline(h=seq(0,1.2*tickz,1.2*tickz/20), lwd=0.5, col="gray" )
            axis(side=2, at=seq(0,1.2*tickz,1.2*tickz/20), labels=yLabel, las=2, cex=0.3, tck=0)
            legend("topleft",c(paste("H1=",nm[1],sep=""),paste("H2=",nm[2],sep=""),paste("H3=",nm[3],sep=""),paste("H4=",nm[4],sep="")), lty = c(1,1,1,1), lwd=c(20,20,20,20), col = c("darkgreen","darkorange1","darkmagenta","firebrick1"),bty="n",ncol=2)

            dev.off()

                                        #Print estimated error rates on a file

            fileOut = paste(dirName,"/errorRates",combNames[ idComb ],".txt",sep="")
            cat(err1,"\n",err2,"\n",err3,"\n",err4,file=fileOut)
        }
        }


        
    }

    

### D WITH NO ERROR CORRECTION AND USING ALL TRANSITIONS
        message("Running Observed D\n",appendLF=FALSE)
        flush.console()
        #fileOut=paste(out,combNames[ idComb ],".Observed",".txt",sep="",collapse="")
        result5 = getJackKnife(outData,ABBAname=ABBA,BABAname=BABA,BBAAname=BBAA)
    if(idComb==1){
        str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBlocks\tH1\tH2\tH3\tH4")
        str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result5$thetaN,result5$thetaJack,result5$varJack,result5$Z,result5$pv,result5$nABBA,result5$nBABA,result5$numBlock,nm[1],nm[2],nm[3],nm[4])
        cat(str,str2,file=FILEOBS,sep="\n")
    }    else{
        str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result5$thetaN,result5$thetaJack,result5$varJack,result5$Z,result5$pv,result5$nABBA,result5$nBABA,result5$numBlock,nm[1],nm[2],nm[3],nm[4])
        cat(str2,file=FILEOBS,sep="\n",append=TRUE)
    }

### D WITH NO ERROR CORRECTION AND REMOVING ANCIENT TRANSITIONS
        message("Running Transition Removal\n",appendLF=FALSE)
        flush.console()

        #fileOut=paste(out,combNames[ idComb ],".TransRem",".txt",sep="",collapse="")

        result6 = getJackKnife(outData,ABBAname=ABBAtr,BABAname=BABAtr,BBAAname=BBAA)
        if(idComb==1){
            str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBlocks\tH1\tH2\tH3\tH4")
            str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result6$thetaN,result6$thetaJack,result6$varJack,result6$Z,result6$pv,result6$nABBA,result6$nBABA,result6$numBlock,nm[1],nm[2],nm[3],nm[4])
            cat(str,str2,file=FILETRANS,sep="\n")
            }  else{
                str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%s\t%s",result6$thetaN,result6$thetaJack,result6$varJack,result6$Z,result6$pv,result6$nABBA,result6$nBABA,result6$numBlock,nm[1],nm[2],nm[3],nm[4])
                cat(str2,file=FILETRANS,sep="\n",append=TRUE)
            }

### output fancy table
    str0=sprintf(" Combination %s %s %s %s\n", nm[1], nm[2], nm[3], nm[4])
    str1=sprintf("  Mode\t\t|Dstat\t\t|sd(Dstat)\t|Djack\t\t|Zscore\t|Pvalue\n")
    str2=sprintf("Observed\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result5$thetaN,sqrt(result5$varJack),result5$thetaJack,result5$Z,result5$pv)
    if(erCor==1){
        str3=sprintf("Err Corr\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result1$thetaN,sqrt(result1$varJack),result1$thetaJack,result1$Z,result1$pv)}
    str4=sprintf("No Trans\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result6$thetaN,sqrt(result6$varJack),result6$thetaJack,result6$Z,result6$pv)
    if(erCor==1){
        str5=sprintf("Err Corr\t|\t\t|\t\t|\t\t|\t|\t\n")
        str6=sprintf("   and\t\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result3$thetaN,sqrt(result3$varJack),result3$thetaJack,result3$Z,result3$pv)
        str7=sprintf("No Trans\t|\t\t|\t\t|\t\t|\t|\t\n")}

    cat("--- Table of Results --- Combination ", idComb," out of ", numComb," ---\n")
    #cat("--- Time Spent ",(proc.time() - ptm2)[3]," sec --- Estimated time left ", (proc.time() - ptm2)[3]*numComb/idComb - (proc.time()-ptm)[3], " sec ---\n"    )
    cat("--- H1=",nm[1]," H2=",nm[2]," H3=",nm[3]," H4=",nm[4]," ---\n")
    cat("---------------------------------------------------------------------------------\n")
    cat(str1)
    cat("---------------------------------------------------------------------------------\n")
    cat(str2)
    if(erCor==1){
        cat("---------------------------------------------------------------------------------\n")
        cat(str3)}
    cat("---------------------------------------------------------------------------------\n")
    cat(str4)
    cat("---------------------------------------------------------------------------------\n")
    if(erCor==1){
        cat(str5)
        cat(str6)
        cat(str7)
        cat("---------------------------------------------------------------------------------\n")}
}
cat("--- Total time ",(proc.time() - ptm2)[3] + finalptm," sec \n"    )
if(addErrors){
    messageEnd = sprintf("plots with effect of removed errors and D statistic files for all the removed errors are in folder %s\n",dirName)
    cat(messageEnd)
}
