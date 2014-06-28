#small functions to calculate tajima based on the wiki page

##for a,b,c,e functions n is the number of chr in the sample
a1 <- function(n) sum(1/(1:(n-1)))
a2 <- function(n) sum(1/(1:(n-1))^2)

b1 <- function(n) (n+1)/(3*(n-1))
b2 <- function(n) (2*(n^2+n+3))/(9*n*(n-1))

c1 <- function(n) b1(n)-1/a1(n)
c2 <- function(n) b2(n)-(n+2)/(a1(n)*n)+a2(n)/a1(n)^2

e1 <- function(n) c1(n)/a1(n)
e2 <- function(n) c2(n)/(a1(n)^2+a2(n))


##assumes the sfs contains the 2 invariable categories. lengths should be 2*n+1,n is number of diploid samples
pair <- function(sfs,n){
    if(length(sfs)!=2*n+1)
        stop("Problem with dims")
    scal <- (0:(2*n))*((2*n):0)
    sum(sfs*scal/choose(2*n,2))
}


pair.fold <- function(sfs,n){
    if(length(sfs)!=2*n+1)
        stop("Problem with dims")
    scal <- (0:(2*n))*((2*n):0)
    sum(sfs*scal/choose(2*n,2))
}




##input is thetaD and number of segregating sites
tajima.d <- function(pairwise,S,n){
    top <- pairwise-S/a1(n)
    bot <- sqrt(e1(n)*S+e2(n)*S*(S-1))
    top/bot
}


##tajima.d(pair(b,20),sum(b),n=40)
