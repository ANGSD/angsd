options(warn=1)
bases<-c("A","C","G","T")
b <- c(bases,"N")

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
args<-list(angsdFile = NULL,
           file1=FALSE,
           file2=FALSE,
           file3=FALSE,
           file4=FALSE,
           out="errorEst",
           erCor=0,
           addErrors=0,
           nIter=100,
           main="",
           maxErr=0.02
           )

#if no argument aree given prints the need arguments and the optional ones with default
des<-list(angsdFile="output angsdFile (no .abbababa2 extension) from 'angsd -doAbbababa2 1' command",
          file1="the ancError File of population H1",
          file2="the ancError File of population H2",
          file3="the ancError File of population H3",
          file4="the ancError File of population H4",
          out="Name of the out files",
          erCor="1 perform the error correction",
          addErrors="Set to 1 to add to the output addition of type-spec errors in the range [-0.005,0.005]",
          nIter="Number of optimazation attemps",
          maxErr="maximum allowed error rate"
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

cat("----------\nangsdFile: ",angsdFile,"  file1:",file1," file2:",file2," file3:",file3," out:",out," addErrors:"," erCor",erCor, addErrors," nIter:",nIter,"\nmaxErr:",maxErr,"\n-----------\n")

b<-c("A","C","G","T","N")

nIter<-as.integer(nIter)
maxErr=as.numeric(maxErr)

readTable <- function(file){
    r <- colSums(as.matrix(read.table(file)))
    r <- matrix(data=r,nrow=1,ncol=length(r))
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
        return(rbind(res,conv$par))
    }
}

buildMat <- function(res){
    cont=1;
    resMat=matrix(nrow=4,ncol=4)
    for(rws in 1:4){
        for(cls in 1:4){
            if(rws == cls){
                resMat[rws,cls] = 1-res[1,(rws-1)*3+1]-res[1,(rws-1)*3+2]-res[1,(rws-1)*3+3];
            }
            else{
                resMat[rws,cls] = res[1,cont]
                cont=cont+1;
            }
        }
    }
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

ABBA<-c("X0110","X0220","X0330","X1001","X1221","X1331","X2002","X2112","X2332","X3003","X3113","X3223")
BABA<-c("X0101","X0202","X0303","X1010","X1212","X1313","X2020","X2121","X2323","X3030","X3131","X3232")
ABBAtr<-c("X0110","X0330","X1001","X1221","X2112","X2332","X3003","X3223")
BABAtr<-c("X0101","X0303","X1010","X1212","X2121","X2323","X3030","X3232")
BBAA<-c("X0011","X0022","X0033","X1100","X1122","X1133","X2200","X2211","X2233","X3300","X3311","X3322")

getJackKnife <- function(outData,finalInv,printData=0,ABBAname,BABAname,BBAAname){
    
    seenSites = sum(as.numeric(outData[,6]))
    weigth = as.numeric(outData[,6])
    
    outData = outData[,-c(1,2,3,4,5,6)]
    Edata <- numeric(prod(dim(outData)))
    dim(Edata)<- dim(outData)
    L = nrow(outData)
    totAbba = 0
    totBaba = 0
    totBbaa = 0
    num = rep(0,L)
    den = rep(0,L)
                                            
    for(blk in 1:L){
        
        x = as.numeric(outData[blk,])
        Edata[blk,] = finalInv %*% x
        x = Edata[blk,]
        names(x) = colnames(outData)
        num[blk] = sum(as.numeric(x[ABBAname])) - sum(as.numeric(x[BABAname]))
        den[blk] = sum(as.numeric(x[ABBAname])) + sum(as.numeric(x[BABAname]))
        totAbba = totAbba + sum(as.numeric(x[ABBAname]))
        totBaba = totBaba + sum(as.numeric(x[BABAname]))
        totBbaa = totBbaa + sum(as.numeric(x[BBAAname]))
    }
    
    #block jack knife estimation
    remIdx = (is.nan(num) | is.nan(den)) | (is.na(num) | is.na(den))
    
    num = num[!remIdx]
    den = den[!remIdx]
    
    weigth = weigth[!remIdx]/sum(weigth[!remIdx]) #block weights
    L = length(num)
    thetaN = sum(num)/sum(den) #D statistics calculated without jackknife
    errorCount = 0
    
    for(i in 1:L)
        if(num[i]>den[i])
            errorCount = errorCount + 1
      
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
    return(list(thetaN=thetaN,thetaJack=thetaJack,varJack=varJack,Z=Z,pv=pv,nABBA=totAbba,nBABA=totBaba,nBBAA=totBbaa))
}

angsdFile = paste(angsdFile,".abbababa2",sep="")
outData <- read.table(angsdFile,header=T,as.is=T,sep="")
erCor=as.numeric(erCor)
if(erCor==1){#ERROR CORRECTED D

    nInd1<-1
    nInd2<-1
    nInd3<-1
    nInd4<-1

    res1=diag(c(1,1,1,1))
    res2=diag(c(1,1,1,1))
    res3=diag(c(1,1,1,1))
    res4=diag(c(1,1,1,1))  

    resMat = list();
    resMat[[1]] = diag(c(1,1,1,1))
    resMat[[2]] = diag(c(1,1,1,1))
    resMat[[3]] = diag(c(1,1,1,1))
    resMat[[4]] = diag(c(1,1,1,1))
    
if(file1!=FALSE){
    r1 = readTable(file1)
    res1 = getFromErrFile(r1,res1,maxErr,nInd1,logLike)
    resMat[[1]] = buildMat(res1)
}
if(file2!=FALSE){
    r2 = readTable(file2)
    res2 = getFromErrFile(r2,res2,maxErr,nInd2,logLike)
    resMat[[2]] = buildMat(res2)
    }
if(file3!=FALSE){
    r3 = readTable(file3)
    res3 = getFromErrFile(r3,res3,maxErr,nInd3,logLike)
    resMat[[3]] = buildMat(res3)
}
if(file4!=FALSE){
    r4 = readTable(file4)
    res4 = getFromErrFile(r4,res4,maxErr,nInd4,logLike)
    resMat[[4]] = buildMat(res4)
}    
    

errMat = getErrMat(resMat)

finalInv = buildInv(errMat)

result1 = getJackKnife(outData,finalInv,printData=1,ABBAname=ABBA,BABAname=BABA,BBAAname=BBAA)

fileOut = paste(out,".ErrorCorr",".txt",sep="")

str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBBAA")
str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",result1$thetaN,result1$thetaJack,result1$varJack,result1$Z,result1$pv,result1$nABBA,result1$nBABA,result1$nBBAA)
cat(str,str2,file=fileOut,sep="\n")

result3 = getJackKnife(outData,finalInv,printData=1,ABBAname=ABBAtr,BABAname=BABAtr,BBAAname=BBAA)

fileOut = paste(out,".TransRemErrorCorr",".txt",sep="")

str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBBAA")
str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",result3$thetaN,result3$thetaJack,result3$varJack,result3$Z,result3$pv,result3$nABBA,result3$nBABA,result3$nBBAA)
cat(str,str2,file=fileOut,sep="\n")

}

### write in files the effect of added error to transition FROM --> TO
### NEEDS TO BE FINISHED
addErr = c(0,0.001,0.002,0.003,0.004,0.005)
LErr = length(addErr)
letters = c("A","C","G","T")

if(addErrors){###########
print("-----------------------------------------------")
print("Add error to transitions")    
for(FROM in 1:4){
    for(TO in 1:4){
        if(FROM!=TO){
            print(sprintf("--changing transitions from %s to %s",letters[FROM],letters[TO]))
            fileOut = paste(out,"Error",letters[FROM],"To",letters[TO],".txt",sep="")
            for(id in 1:3){
                for(er in 1:LErr){
                    newMat = resMat
                    newMat[[id]][FROM,TO] = newMat[[id]][FROM,TO] + addErr[er]
                    newMat[[id]][FROM,FROM] = newMat[[id]][FROM,FROM] - addErr[er]
                    newErrMat = getErrMat(newMat)
                    newInv = buildInv(newErrMat)
                    result = getJackKnife(outData,newInv,printData=0,ABBAname=ABBA,BABAname=BABA,BBAAname=BBAA)
                    
                    if(file.exists(fileOut))
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv,result2$nABBA,result2$nBABA,result2$nBBAA),fileOut,sep="\t",append=T,ncolumns=10)
                    if(!file.exists(fileOut)){
                        file.create(fileOut, showWarnings = FALSE)
                        str = sprintf("NumInd\taddErr\tD\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBBAA")
                        write(str,fileOut,append=T)
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv,result$nABBA,result$nBABA,result$nBBAA),fileOut,sep="\t",append=T,ncolumns=10)
                    }
                }    
            }
        }
    }
}
        
for(FROM in 1:4){
    for(TO in 1:4){
        if(FROM!=TO){
            print(sprintf("--changing transitions from %s to %s",letters[FROM],letters[TO]))
            fileOut = paste(out,"NoTransError",letters[FROM],"To",letters[TO],".txt",sep="")
            for(id in 1:3){
                for(er in 1:LErr){
                    newMat = resMat
                    newMat[[id]][FROM,TO] = newMat[[id]][FROM,TO] + addErr[er]
                    newMat[[id]][FROM,FROM] = newMat[[id]][FROM,FROM] - addErr[er]
                    newErrMat = getErrMat(newMat)
                    newInv = buildInv(newErrMat)
                    result = getJackKnife(outData,newInv,printData=0,ABBAname=ABBAtr,BABAname=BABAtr,BBAAname=BBAA)
                    
                    if(file.exists(fileOut))
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv,result$nABBA,result$nBABA,result$nBBAA),fileOut,sep="\t",append=F,ncolumns=10)
                    if(!file.exists(fileOut)){
                        file.create(fileOut, showWarnings = FALSE)
                        str = sprintf("NumInd\taddErr\tD\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBBAA")
                        write(str,fileOut,append=T)
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv,result$nABBA,result$nBABA,result$nBBAA),fileOut,sep="\t",append=T,ncolumns=10)
                   }
                }
            }
        }
    }     
}
}#####end if(FALSE)

### D WITH NO ERROR CORRECTION AND USING ALL TRANSITIONS
fileOut=paste(out,".Observed",".txt",sep="",collapse="")

outData <- read.table(angsdFile,header=T,as.is=T,sep="")
result5 = getJackKnife(outData,diag(rep(1,256)),printData=1,ABBAname=ABBA,BABAname=BABA,BBAAname=BBAA)

str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBBAA")
str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",result5$thetaN,result5$thetaJack,result5$varJack,result5$Z,result5$pv,result5$nABBA,result5$nBABA,result5$nBBAA)
cat(str,str2,file=fileOut,sep="\n")

### D WITH NO ERROR CORRECTION AND REMOVING ANCIENT TRANSITIONS
fileOut=paste(out,".RemTrans",".txt",sep="",collapse="")


outData <- read.table(angsdFile,header=T,as.is=T,sep="")
result6 = getJackKnife(outData,diag(rep(1,256)),printData=1,ABBAname=ABBAtr,BABAname=BABAtr,BBAAname=BBAA)

str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue\tnABBA\tnBABA\tnBBAA")
str2 = sprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",result6$thetaN,result6$thetaJack,result6$varJack,result6$Z,result6$pv,result6$nABBA,result6$nBABA,result6$nBBAA)
cat(str,str2,file=fileOut,sep="\n")

### output fancy table
str1=sprintf("  Mode\t\t|Dstat\t\t|sd(Dstat)\t|Djack\t\t|Zscore\t|Pvalue\n")
str2=sprintf("Observed\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result5$thetaN,sqrt(result5$varJack),result5$thetaJack,result5$Z,result5$pv)
if(erCor==1){
str3=sprintf("Err Corr\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result1$thetaN,sqrt(result1$varJack),result1$thetaJack,result1$Z,result1$pv)}
str4=sprintf("No Trans\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result6$thetaN,sqrt(result6$varJack),result6$thetaJack,result6$Z,result6$pv)
if(erCor==1){
str5=sprintf("Err Corr\t|\t\t|\t\t|\t\t|\t|\t\n")
str6=sprintf("   and\t\t|%.3e\t|%.3e\t|%.3e\t|%.3f\t|%.1e\n",result3$thetaN,sqrt(result3$varJack),result3$thetaJack,result3$Z,result3$pv)
str7=sprintf("No Trans\t|\t\t|\t\t|\t\t|\t|\t\n")}

cat("--- Table of Results ---\n")
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
