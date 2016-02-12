
options(warn=1)
bases<-c("A","C","G","T")
b <- c(bases,"N")
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
  cat("->  Needed arguments:\n")
  mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("->  Optional arguments (defaults):\n")
  mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}

## choose your parameters and defaults
## NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(angsdFile = NULL,
           file1="NoFile1",
           file2="NoFile2",
           file3="NoFile3",
           out="errorEst",
           erCor=0,
           addErrors=0,
           nIter=100,
           main="Error rate using an outgroup and a high quality genome",
           maxErr=0.02
           )
#if no argument aree given prints the need arguments and the optional ones with default
des<-list(angsdFile="output angsdFile (w extension) from 'angsd -doAbbababa2 1' command",
          file1="the ancError File of population H1",
          file2="the ancError File of population H2",
          file3="the ancError File of population H3",
          out="Name of the out files",
          erCor="1 perform the error correction",
          addErrors="Set to 1 to add to the output addition of type-spec errors in the range [0,0.005]",
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
#angsdFile="errortest01WG";doError=0;file1="/ricco/genomes/humanAnalysis/errorRates/output/Clovis.flt.sort.rmdup.realign.mdQ20q30V1.ancError";file2="/ricco/genomes/humanAnalysis/errorRates/output/HGDP00877.hg19.flt.sort.rmdup.realign.mdQ20q30V1.ancError";file3="/ricco/genomes/humanAnalysis/errorRates/output/HGDP00521.hg19.flt.sort.rmdup.realign.mdQ20q30V1.ancError"; out="errortest01";nIter=100;main="dfsf";maxErr=0.02;



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
    decomp = svd(res);
    return(decomp$v %*% diag(1/decomp$d) %*% t(decomp$u));       
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

getJackKnife <- function(outData,finalInv,printData=0,ABBAname,BABAname){
    
    seenSites = sum(as.numeric(outData[,6]))
    weigth = as.numeric(outData[,6])
    
    outData = outData[,-c(1,2,3,4,5,6)]
    Edata <- numeric(prod(dim(outData)))
    dim(Edata)<- dim(outData)
    L = nrow(outData)
    totAbba = rep(0,L)
    totBaba = rep(0,L)
    num = rep(0,L)
    den = rep(0,L)
                                            
    for(blk in 1:L){
        
        x = as.numeric(outData[blk,])
        Edata[blk,] = finalInv %*% x
        x = Edata[blk,]
        names(x) = colnames(outData)
        num[blk] = sum(as.numeric(x[ABBAname])) - sum(as.numeric(x[BABAname]))
        den[blk] = sum(as.numeric(x[ABBAname])) + sum(as.numeric(x[BABAname]))
        
    }
    
    #block jack knife estimation
    remIdx = (is.nan(num) | is.nan(den)) | (is.na(num) | is.na(den))
    
    num = num[!remIdx]
    den = den[!remIdx]
    
    weigth = weigth[!remIdx]/sum(weigth[!remIdx])
    L = length(num)
    thetaN = sum(num)/sum(den)
    errorCount = 0
    
    for(i in 1:L)
        if(num[i]>den[i])
            errorCount = errorCount + 1
    
    
    thetaJStar <- rep(0,L) #partial estimators
    thetaJack <- rep(0,L)  #Block Jack knife estimator
    pseudo <- 1/weigth
    thetaJStar <- sapply( 1:L, function(x) sum(as.numeric(num[-x]))/sum(as.numeric(den[-x])) )
    thetaJack <- L*thetaN - sum((1-weigth) * thetaJStar)        
    meanJack <- mean(thetaJStar)
    thetaTilde <- pseudo*thetaN-(pseudo-1) * thetaJStar
    varJack <- 1/L * sum( 1/(pseudo-1) * (thetaTilde - thetaJack)^2 )       
    Z <- thetaN / sqrt(varJack)
    pv = 2*min(pnorm(Z,0,1),1-pnorm(Z,0,1))
    if(printData){
        print(c("thetaN",thetaN))
        print(c("theta JF",thetaJack))
        print(c("SD D_Jack",sqrt(varJack)))
        print(c("Z",Z))
        print(c("p-value",pv))
        print(c("visited sites",seenSites))
        if(errorCount > 0)
            print(sprintf("Warning: you got %d times that Num>Den. This can  happen when you apply error correction on alleles combination with low probability. If this happen too many times, may be error correction is unnecessary.",errorCount))
    }
    return(list(thetaN=thetaN,thetaJack=thetaJack,varJack=varJack,Z=Z,pv=pv))
}





angsdFile = paste(angsdFile,".abbababa2",sep="")
outData <- read.table(angsdFile,header=T,as.is=T,sep="")

if(erCor){
print("-----------------------------------------------")    
print("Estimation with error correction and no Ancient Transition removal")

r1 = readTable(file1)
r2 = readTable(file2)
r3 = readTable(file3)

nInd1<-1
nInd2<-1
nInd3<-1

res1 = NULL
res2 = NULL
res3 = NULL
res1 = getFromErrFile(r1,res1,maxErr,nInd1,logLike)
res2 = getFromErrFile(r2,res2,maxErr,nInd2,logLike)
res3 = getFromErrFile(r3,res3,maxErr,nInd3,logLike)

resMat = list();

resMat[[1]] = buildMat(res1)
resMat[[2]] = buildMat(res2)
resMat[[3]] = buildMat(res3)
resMat[[4]] = diag(rep(1,4))

errMat = getErrMat(resMat)

finalInv = buildInv(errMat)

result1 = getJackKnife(outData,finalInv,printData=1,ABBAname=ABBA,BABAname=BABA)

fileOut = paste(out,"ErrorCorr",".txt",sep="")

#if(file.exists(fileOut))
#    write(c(result1$thetaN,result1$thetaJack,result1$varJack,result1$Z,result1$pv),fileOut,sep="\t",append=T)
#if(!file.exists(fileOut)){
#    file.create(fileOut, showWarnings = FALSE)
#    str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue")
#    write(str,fileOut,append=T)
#    write(c(result1$thetaN,result1$thetaJack,result1$varJack,result1$Z,result1$pv),fileOut,sep="\t",append=T)
#}

str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue")
str2 = sprintf("%f\t%f\t%f\t%f\t%f",result1$thetaN,result1$thetaJack,result1$varJack,result1$Z,result1$pv)
cat(str,str2,file=fileOut,sep="\n")




print("-----------------------------------------------")
print("Error correction and removal of ancient transition")
result2 = getJackKnife(outData,finalInv,printData=1,ABBAname=ABBAtr,BABAname=BABAtr)

fileOut = paste(out,"ErrorCorrNoTrans",".txt",sep="")
#if(!file.exists(fileOut)){
#    file.create(fileOut, showWarnings = FALSE)
#    str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue")
#    write(str,fileOut,append=T)
#    write(c(result2$thetaN,result2$thetaJack,result2$varJack,result2$Z,result2$pv),fileOut,sep="\t",append=T,ncolumns=5)
#}
#if(file.exists(fileOut))
#    write(c(result2$thetaN,result2$thetaJack,result2$varJack,result2$Z,result2$pv),fileOut,sep="\t",append=T,ncolumns=5)

str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue")
str2 = sprintf("%f\t%f\t%f\t%f\t%f",result2$thetaN,result2$thetaJack,result2$varJack,result2$Z,result2$pv)
cat(str,str2,file=fileOut,sep="\n")


}






### write in files the effect of added error to transition FROM --> TO

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
                    result = getJackKnife(outData,newInv,printData=0,ABBAname=ABBA,BABAname=BABA)
                    
                    if(file.exists(fileOut))
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv),fileOut,sep="\t",append=T,ncolumns=7)
                    if(!file.exists(fileOut)){
                        file.create(fileOut, showWarnings = FALSE)
                        str = sprintf("NumInd\taddErr\tD\tJK-D\tV(JK-D)\tZ\tpvalue")
                        write(str,fileOut,append=T)
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv),fileOut,sep="\t",append=T,ncolumns=7)
                    }

#str = sprintf("NumInd\taddErr\tmean(D)\tJK-D\tV(JK-D)\tZ\tpvalue")
#str2 = sprintf("%d\t%f\t%f\t%f\t%f\t%f\t%f",id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv)
#cat(str,str2,file=fileOut,sep="\n")
                    
                }    
            }
        }
    }
}


print("-----------------------------------------------")
print("Remove ancient transitions G->A + C->T  AND  Add error to transitions")
        
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
                    result = getJackKnife(outData,newInv,printData=0,ABBAname=ABBAtr,BABAname=BABAtr)
                    
                    if(file.exists(fileOut))
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv),fileOut,sep="\t",append=F,ncolumns=7)
                    if(!file.exists(fileOut)){
                        file.create(fileOut, showWarnings = FALSE)
                        str = sprintf("NumInd\taddErr\tD\tJK-D\tV(JK-D)\tZ\tpvalue")
                        write(str,fileOut,append=T)
                        write(c(id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv),fileOut,sep="\t",append=T,ncolumns=7)
                   }


#str = sprintf("NumInd\taddErr\tmean(D)\tJK-D\tV(JK-D)\tZ\tpvalue")
#str2 = sprintf("%d\t%f\t%f\t%f\t%f\t%f\t%f",id,addErr[er],result$thetaN,result$thetaJack,result$varJack,result$Z,result$pv)
#cat(str,str2,file=fileOut,sep="\n")
                    
                    
                }
            }
        }
    }     
}
}#####end if(FALSE)

print("-----------------------------------------------")
print("D-statistic calculated without Error correction")

fileOut=paste(out,"Std",".txt",sep="",collapse="")

outData <- read.table(angsdFile,header=T,as.is=T,sep="")
result5 = getJackKnife(outData,diag(rep(1,256)),printData=1,ABBAname=ABBA,BABAname=BABA)

#if(file.exists(fileOut))
#    write(c(result5$thetaN,result5$thetaJack,result5$varJack,result5$Z,result5$pv),fileOut,sep="\t",append=T)
#if(!file.exists(fileOut)){
#    file.create(fileOut, showWarnings = FALSE)
#    str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue")
#    write(str,fileOut,append=T)
#    write(c(result5$thetaN,result5$thetaJack,result5$varJack,result5$Z,result5$pv),fileOut,sep="\t",append=T)
#}


str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue")
str2 = sprintf("%f\t%f\t%f\t%f\t%f",result5$thetaN,result5$thetaJack,result5$varJack,result5$Z,result5$pv)
cat(str,str2,file=fileOut,sep="\n")

print("-----------------------------------------------")
print("D-statistic calculated w/o Error correction and transitions")

fileOut=paste(out,"NoErrorNoTrans",".txt",sep="",collapse="")

outData <- read.table(angsdFile,header=T,as.is=T,sep="")
result6 = getJackKnife(outData,diag(rep(1,256)),printData=1,ABBAname=ABBAtr,BABAname=BABAtr)

#if(file.exists(fileOut))
#    write(c(result6$thetaN,result6$thetaJack,result6$varJack,result6$Z,result6$pv),fileOut,sep="\t",append=T)
#if(!file.exists(fileOut)){
#    file.create(fileOut, showWarnings = FALSE)
#    str = sprintf("D\tJK-D\tV(JK-D)\tZ\tpvalue")
#    write(str,fileOut,append=T)
#    write(c(result6$thetaN,result6$thetaJack,result6$varJack,result6$Z,result6$pv),fileOut,sep="\t",append=T)
#}

str = sprintf("mean(D)\tJK-D\tV(JK-D)\tZ\tpvalue")
str2 = sprintf("%f\t%f\t%f\t%f\t%f",result6$thetaN,result6$thetaJack,result6$varJack,result6$Z,result6$pv)
cat(str,str2,file=fileOut,sep="\n")
