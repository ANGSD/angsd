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

funner <- function(from,to,st.fix,ref.fix,posi=T){
  r <- c()
  for(i in from:to)
      if(posi)
          r <- rbind(r,norm(colSums(cnts[p1==i&ref==ref.fix&st==st.fix,])))
      else
          r <- rbind(r,norm(colSums(cnts[p2==i&ref==ref.fix&st==st.fix,])))
  colnames(r) <- c("A","C","G","T") 
  r
}

pdf("test.pdf",width=28,height=28)
par(mfrow=c(4,4))
for(r in 0:3)
    for(s in 0:1)
        for(p in c(T,F)){
            res <- funner(from=0,to=20,st.fix=s,ref.fix=r,posi=p)*100
            matplot(res[,-(r+1)],lwd=c(2),type='b',lty=1,pch=c("A","C","G","T")[-(r+1)],main=paste("ref:",bases[r+1],"strand: ",s," p:",p))
        }

dev.off()
