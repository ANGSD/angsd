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

funner <- function(from,to,st.fix,ref.fix){
  r <- c()
  for(i in from:to)
    r <- rbind(r,norm(colSums(cnts[p1==i&ref==ref.fix&st==st.fix,])))
  
  r
}

funner(from=0,to=10,st.fix=0,ref.fix=0)

        
