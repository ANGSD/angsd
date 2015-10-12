pplot<-function(x,ylab="YRI",xlab="CEU",pal,...){
    x[x<1e-9]<-1e-9
    p2<-x
    vek<-sort(c(1,0.05,0.01,0.005,0.001,0.0005,10^-(4:10),5*10^-(4:10)),decreasing=T)
    ##vek<-vek[1:which.max(min(x,na.rm=T)>vek)]
    vek<-vek[(which.max(!(max(x,na.rm=T)<vek))-1):which.max(min(x,na.rm=T)>vek)]


    v<-length(vek)
    for(tal in 1:v)
        p2[x<vek[tal]]<-v-tal


    if(missing(pal))
        ccol<-heat.colors(v+2)
    else
        ccol<-pal(v+2)
    lattice::levelplot(p2,region=T,colorkey=list(at=1:v+1,label=list(at=1:v+1,lab=as.character(rev(vek)))),col.regions=ccol,cut=v,ylab=ylab,xlab=xlab,scales= list(x = list(rot = 45)),...)
}

color.palette <- function(steps, n.steps.between=NULL, ...){

    if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
    if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")

    fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[,fill.steps] <- col2rgb(steps)

    for(i in which(n.steps.between>0)){
        col.start=RGB[,fill.steps[i]]
        col.end=RGB[,fill.steps[i+1]]
        for(j in seq(3)){
            vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]
            RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
        }
    }

    new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}
if(FALSE){

    require(lattice)
    source("plot2dSFS.R")
    
    pal <- color.palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","red4","darkred","black"), space="rgb")

  est<-as.matrix(exp(read.table("/ricco/home/home/thorfinn/sfs_greenland/ext_accessible2/chb.gr.2dsfs.ml")))
    rownames(est)<-1:nrow(est)-1 # Greenland
    colnames(est)<-1:ncol(est)-1 # China

    minVal<-1e-9
    est[est<minVal]<-minVal
est0<-est
    est0[1,1]<-NA
    pdf("~/public/albrecht/open/tmp.pdf")
    pplot(est,pal=pal,ylab="CHB",xlab="Greenland",main="2D site frequency spectrum")
    pplot(est0,pal=pal,ylab="CHB",xlab="Greenland",main="2D site frequency spectrum")
    dev.off()


    ## come to dadi
    pol<-"folded"
    sfs<-est*1e7
    n1=nrow(sfs);
    n2=ncol(sfs);
    ns=paste(n1,n2,collapse=' ');
    ns=paste(ns,pol,collapse=' ');
    ## Convert 2D sfs to dadi array format
    dadi=NULL;
    for(i in 1:n1){
        dadi=c(dadi,as.numeric(sfs[i,]))
    }
    file<-"tmpDadi"
# Write out dadi format to file with same name as 2D sfs file with .fs appeded
write.table(ns,file=paste(file,'.fs',sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE);
write.table(paste(dadi,collapse=' '),file=paste(file,'.fs',sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);


}
