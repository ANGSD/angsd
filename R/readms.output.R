
#read ms output
#
# The argument is either a file name or a vector of character strings, one 
# string for each line of the output of ms.
#  The function returns a list with some of the following components: 
#       segsites,  times, positions, gam, probs, nsam, nreps
# Example usage:
#  system("./ms 5 2 -s 5 >ms.out")
#  msout <- read.ms.output(file="ms.out")
#   Then for example, msout$gam[[1]] is a haplotype array for the first sample
#  or 
# msout.txt <- system("./ms 5 2 -s 5 -L", intern=TRUE)
#  msout <- read.ms.output(msout.txt)
# msout$time[1,1] is then the tmrca of the first sample
# msout$time[1,2] is the total tree length of the first sample
# or 
#  msout <- read.ms.output( system("./ms 5 3 -t 3.0 -s 4 -L",intern=TRUE))
# This function is derived from code first written by Dan Davison.

read.ms.output <- function( txt=NA, file.ms.output=NA ) {
    
    if( !is.na(file.ms.output) ) txt <- scan(file=file.ms.output,
       what=character(0), sep="\n", quiet=TRUE)
    if( is.na(txt) ){
    	print("Usage: read.ms.output(txt), or read.ms.output(file=filename)")
    	return()
    	}
    nsam <- as.integer( strsplit(txt[1], split=" ")[[1]][2] )
    ndraws <- as.integer( strsplit( txt[1], split=" ")[[1]][3] )

    h <- numeric()
    result <- list()
    gamlist <- list()
    positions <- list()

    marker <- grep("prob",txt)
    probs <- sapply(strsplit(txt[marker], split=":"), function(vec) as.numeric(vec[2]))
    marker <- grep("time",txt)
    times <- sapply(strsplit(txt[marker], split="\t"), function(vec){ as.numeric(vec[2:3])} )

    
    ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
    marker <- grep("segsites", txt)
    stopifnot(length(marker) == ndraws)

    ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
    segsites <- sapply(strsplit(txt[marker], split=" "), function(vec) as.integer(vec[2]) )
    for(draw in seq(along=marker)) {
        if(!(draw %% 100)) cat(draw, " ")
        if(segsites[draw] > 0) {
        	  tpos <- strsplit(txt[marker[draw]+1], split=" ")
        	  positions[[draw]] <- as.numeric( tpos[[1]][ 2:(segsites[draw]+1) ] ) 
            haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
            haplotypes <- strsplit(haplotypes, split="")
            h <- sapply(haplotypes, function(el) c(as.integer(el)))
            ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
            if(segsites[draw] == 1) h <- as.matrix(h)
            ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
            else h <- t(h)
        }
        else {
        	h <- matrix(nrow=nsam, ncol=0)
        	positions[[draw]]<- NA
        }
		 gamlist[[draw]] <- h
        stopifnot(all(dim(h) == c(nsam, segsites[draw]))) 
    }
	cat("\n")
    list(segsites=segsites, gam=gamlist, probs=probs, times=t(times), positions=positions, nsam=nsam, nreps=ndraws ) 
}
