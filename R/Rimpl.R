
readGeno <- function(file=NULL,nind=10,nsites=10){
  ff <- gzfile(file,"rb")
  m<-matrix(readBin(ff,double(),10*nind*nsites),ncol=10*nind,byrow=TRUE)
  close(ff)
  return(m)
}

readBjoint <- function(file=NULL,nind=10,nsites=10){
  ff <- gzfile(file,"rb")
  m<-matrix(readBin(ff,double(),(2*nind+1)*nsites),ncol=(2*nind+1),byrow=TRUE)
  close(ff)
  return(m)
}

##we need a majorminor table
majorminor <- rbind(c(0,1,2,3),c(1,4,5,6),c(2,5,7,8),c(3,6,8,9))+1


#likes <- readGeno(file="pops1.glf.gz",nind=5,nsites=100000) 



##est MAF for single site given the major minor
emFrequency1 <- function(likes1,majmin=c(1,2),start=0.001,iter=10){
  major=majmin[1]
  minor=majmin[2]
  if(class(likes1)!="matrix")
    likes1 <- matrix(likes1,nrow=10)
 
  p <- start
  sum <- c()
  for(i in 1:iter){
    tmp <- c()
    tmp <- rbind(tmp,exp(likes1[majorminor[major,major],])*(1-p)^2)
    tmp <- rbind(tmp,exp(likes1[majorminor[major,minor],])*2*p*(1-p))
    tmp <- rbind(tmp,exp(likes1[majorminor[minor,minor],])*p^2)
 #   print((tmp))
    ws <- sum((tmp[2,]+2*tmp[3,])/(2*colSums(tmp)))
  #  cat("ws:",ws,"\n")
    p <- ws/ncol(likes1)
   # print(p)
   # stop("asdf")
  }
  return (p)
}

#emFrequency1(likes1=likes[1,],iter=10)

estMajorMinor1 <- function(likes1){
  likes1 <- matrix(likes1,nrow=10)

  lmax <- -10000000
  totallike <- c()
  choiceMajor <- c()
  choiceMinor <- c()
  for(imajor in 1:3)
    for(iminor in (imajor+1):4){
     ## cat(imajor,iminor,"\n")
      tmp <- c()
      tmp <- rbind(tmp,log(0.25)+likes1[majorminor[imajor,imajor],])
      tmp <-rbind(tmp, log(0.5)+likes1[majorminor[imajor,iminor],])
      tmp <- rbind(tmp,log(0.25)+likes1[majorminor[iminor,iminor],])
      totallike <- sum(apply(tmp,2,function(x) log(sum(exp(x)))))
      #print(tmp)
      #cat(totallike,"\n")
      if(totallike>lmax){
        lmax <- totallike
        choiceMajor <- imajor
        choiceMinor <- iminor
      }
      
    }
  majmin <- c(choiceMajor,choiceMinor)
  MAF <- emFrequency1(likes1=likes1,majmin=majmin)
  if(MAF>0.5)
    return(list(MAJMIN=rev(majmin),maf=1-MAF))
  else
    return(list(MAJMIN=majmin,maf=MAF))
}


##estMajorMinor1(likes[1,])

#majminmaf <- apply(likes,1,estMajorMinor1)
#head(t(matrix(unlist(majminmaf),nrow=3)))

estsfs1 <- function(likes1,anc){
  nind <- length(likes1)/10 ##number of inds
  ressfs <- rep(0,nind*2+1) ##result array, sum of 3 derived in regular space
  for(i in 1:4){
    tmpsfs <- rep(0,nind*2+1)#tmp array for one derived
   
    if(i==anc)
      next #skip if derived equals ancestral
    offs <- c(majorminor[anc,anc],majorminor[anc,i],majorminor[i,i])#the offsets
    totmax <- 0 #underflow stuff maybe not needed
    for(i in 1:nind){#loop through samples
      probs <- likes1[(i-1)*10+offs] #extract the 3 likes of interest
      probs[2] <- probs[2]+log(2) #multiply the hetero with 2
      mymax <- max(probs) #rescale for avoiding small values
      probs <- probs-mymax
      totmax <- totmax+mymax
      probs <- exp(probs)
      if(i==1){#hook for first first sample
        tmpsfs[1:3] <- probs
        next
      }
      for(j in (2*i +1):(3)){
        ##tmpsfs[j] <- tmpsfs[j-2]*probs[3]+tmpsfs[j-1]*probs[2]+tmpsfs[j]*probs[1]
        tmpsfs[j] <- sum(probs*tmpsfs[j:(j-2)]) #same as above
      }
      tmpsfs[2] <- tmpsfs[2] * probs[1] + tmpsfs[1]*probs[2]
      tmpsfs[1] <- tmpsfs[1] *probs[1]
    }
    ressfs <-ressfs+ exp(log(tmpsfs)-lchoose(2*nind,0:(2*nind))+totmax)
  }
  ressfs <- log(ressfs)
  ressfs <- ressfs-max(ressfs)
  return(ressfs)
}
#estsfs1(likes[1,],anc=1)
 
## validate like ./angsd.g++ -sim1 pops1.glf.gz -nInd 5  -outfiles testout -doMaf 2 -realSFS 1
dd <- readBjoint(file="testout.sfs",nind=5,nsites=10000)
likes <- readGeno(file="pops1.glf.gz",nind=5,nsites=10000)
aa<-estsfs1(likes[1,],anc=1)
dd[1,]
aa
