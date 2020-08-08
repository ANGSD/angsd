
#' Calculate Pi, via the sum of site heterozygosities 
#' 
#' @param x a gamete matrix
#' @param replace calculate theta pi based on sampling with replacement formula
#' @export
theta_pi <- function(x, replace=FALSE) {
  k <- colSums(x)
  n <- nrow(x)
  if (replace) {
    p <- k/n
    return(sum(2*p*(1-p)))
  } else {
    return(sum(2*(k/n)*((n-k)/(n-1))))
  }
}

theta_pi2 <- function(k,n, replace=FALSE) {
    if (replace) {
    p <- k/n
    return(sum(2*p*(1-p)))
  } else {
    return(sum(2*(k/n)*((n-k)/(n-1))))
  }
}

#' Calculate Watterson's Theta
#' 
#' @param s number of segregating sites
#' @param n number of samples
#' @export
theta_W <- function(s, n) s/(sum(1/(1:(n-1))))

tajd_num <- function(pi, s, n) {
  a1 <- sum(1/(1:(n-1)))
  pi-s/a1
} 

tajd_denom <- function(s, n) {
  a1 <- sum(1/(1:(n-1)))
  a2 <- sum(1/(1:(n-1))^2)
  b1 <- (n+1)/(3*(n-1))
  b2 <- (2*(n^2 + n + 3))/(9*n*(n-1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n+2)/(a1*n) + a2 / a1^2
  e1 <- c1/a1 
  e2 <- c2/(a1^2 + a2)
  sqrt(e1*s + e2*s*(s-1))
}


#' Calculate Tajima's D
#' 
#' @param pi theta pi, or pairwise differences
#' @param s numbr of segregating sites
#' @param n number of samples
#' @export
tajd <- function(pi, s, n) {
  a1 <- sum(1/(1:(n-1)))
  a2 <- sum(1/(1:(n-1))^2)
  b1 <- (n+1)/(3*(n-1))
  b2 <- (2*(n^2 + n + 3))/(9*n*(n-1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n+2)/(a1*n) + a2 / a1^2
  e1 <- c1/a1 
  e2 <- c2/(a1^2 + a2)
  D <- (pi-s/a1)/sqrt(e1*s + e2*s*(s-1))
  D
}


#' Mutate a tibble of results from parse_ms(), adding summary statistics
#' 
#' @param x tibble of results from parse_ms()
#' @param .n number of samples
#'
#' @export
sample_stats <- function(x, .n) {
  x %>% mutate(theta_pi=map_dbl(gametes, theta_pi), theta_W=theta_W(segsites, .n), D=tajd(theta_pi, segsites, .n))
# D_num=tajD_num(theta_pi, segsites, .n), D_denom=tajD_denom(segsites, .n))
}



source("../../../R/readms.output.R")
a <- read.ms.output(file='../input/msout')
norm<-function(x) x/sum(x)
sfs <- as.numeric(table(unlist(sapply(a$gam,colSums))))

theta_pi2(as.numeric(unlist(sapply(a$gam,colSums))),40)
theta_W(sum(sapply(a$gam,ncol)),40)
tajd(theta_pi2(as.numeric(unlist(sapply(a$gam,colSums))),40),264256,40)


> theta_pi2(as.numeric(unlist(sapply(a$gam,colSums))),40)
[1] 55305.51
> theta_W(sum(sapply(a$gam,ncol)),40)
[1] 62126.09
> tajd(theta_pi2(as.numeric(unlist(sapply(a$gam,colSums))),40),264256,40)
[1] -0.4137441
> 
[1]+  Stopped                 R -q  (wd: ~/angsd/test/haploid_sfs/R)
(wd now: ~/angsd/test/haploid_sfs)
thorfinns-mbp:haploid_sfs fvr124$ cat haploid_sfs/output/
thorfinns-mbp:haploid_sfs fvr124$ cat output/fold.thetaFromSaf.thetas.idx.pestPG 
#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)	Chr	WinCenter	tW	tP	tF	tH	tL	Tajima	fuf	fud	fayh	zeng	nSites
(0,264255)(1,264256)(0,264256)	1	132128	62125.858577	55303.560820	0.000000	0.000000	0.000000	-0.413850	1.676630	2.588579	0.714300	-0.813258	264255
thorfinns-mbp:haploid_sfs fvr124$ cat output/norm.thetaFromSaf.thetas.idx.pestPG 
#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)	Chr	WinCenter	tW	tP	tF	tH	tL	Tajima	fuf	fud	fayh	zeng	nSites
(0,264255)(1,264256)(0,264256)	1	132128	62125.858577	55303.579166	69276.207690	98870.688056	77087.134378	-0.413849	-0.423606	-0.297931	-0.562712	0.195850	264255
thorfinns-mbp:haploid_sfs fvr124$ 
