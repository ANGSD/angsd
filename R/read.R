read3gl <- function(name="../angsdput",nSites=1e6){
  p <- paste0(name,".glf.pos.gz")
  g <- paste0(name,".glf.gz")
  con=gzfile(g,"rb");
  glf <- readBin(con,what="double",3*nSites);
  close(con)
  nSites <- length(glf)/3
  
  con=gzfile(p,"rb")
 
 

  ref <- rep(NA,nSites)
  pos <- rep(NA,nSites)
  maj <- rep(NA,nSites)
  min <- rep(NA,nSites)
  for(i in 1:nSites){
    ref[i] <- readBin(con,what=character(),1)
    pos[i] <- readBin(con,what="integer",1)
    maj[i] <- readChar(con,1)
    min[i] <- readChar(con,1)
  }
  toInt <- function(x)  {
    x[x==""] <- 0
    x[x=="\001"] <- 1
    x[x=="\002"] <- 2
    x[x=="\003"] <- 3
    as.integer(x)
  }
  maj <- toInt(maj)
  min <- toInt(min)
  close(con)
  return(list(glf=glf,pos=pos+1,ref=ref,maj=maj,min=min))
}

