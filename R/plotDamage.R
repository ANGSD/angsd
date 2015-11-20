bases <- c("A","C","G","T")

norm <- function(x) x/sum(x)
a <- read.table("angsdput.mismatch.gz",he=T)
a <- a[rowSums(a[,-(1:5)])>0,]

p1 <- a[,1]
p2 <- a[,2]
qs <- a[,3]
st <- a[,4]
ref <- a[,5]
cnts <- a[,1+5:8]

#c->t, G->A

funner <- function(st.fix,ref.fix,posi=T){
    r <- c()
    if(posi){
        pos <- p1
    }else{
        pos <- p2
    }
    for(i in min(pos):max(pos))
        r <- rbind(r,norm(colSums(cnts[pos==i&ref==ref.fix&st==st.fix,])))
    
    colnames(r) <- c("A","C","G","T") 
    r
}

pdf("test.pdf",width=28,height=28)
par(mfrow=c(4,4))
for(r in 0:3)
    for(s in 0:1)
        for(p in c(T,F)){
            res <- funner(from=0,to=93,st.fix=s,ref.fix=r,posi=p)*100
            matplot(res[,-(r+1)],lwd=c(2),type='b',lty=1,pch=c("A","C","G","T")[-(r+1)],main=paste("ref:",bases[r+1],"strand: ",s," p:",p))
        }

dev.off()

##flip

cnts[st==1,4:1] <- cnts[st==1,1:4]



pdf("test2.pdf",width=28,height=28)
par(mfrow=c(4,4))
for(r in 0:3)
    for(s in 0:1)
        for(p in c(T,F)){
            res <- funner(from=0,to=93,st.fix=s,ref.fix=r,posi=p)*100
            matplot(res[,-(r+1)],lwd=c(2),type='b',lty=1,pch=c("A","C","G","T")[-(r+1)],main=paste("ref:",bases[r+1],"strand: ",s," p:",p))
        }

dev.off()

r1<-funner(st.fix=0,ref.fix=0,posi=T)
r2<-funner(st.fix=1,ref.fix=0,posi=F)


r11 <- r1[,-1]
r21 <- r2[,4:1]
r21[,4] <- 1- r21[,4]
r21 <- r21[,-1]

res <- cbind(r11,r21[nrow(r21):1,])
id <- c(paste0(bases[-1],'+'),paste0(bases[-1],'-'))
matplot(res,lwd=c(2),type='b',pch=bases[-1],col=1:6,lty=rep(1:2,each=3))
legend("center",id,fill=1:6)
