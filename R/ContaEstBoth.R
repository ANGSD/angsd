#################
#parse arguments
args<-commandArgs(TRUE)

#names of recognized arguments
PossibleArgNames<-c("counts", "freqs", "outfile", "maxsites", "nthr", "oneCns")

USAGE<-("counts\tSTR\tangsd/misc/contamination -k 1 output. Required.
	freqs\tSTR\tHapMap file. Should be same as angsd/misc/contamination -h. Required.
	outfile\tSTR\tOutput file name. Required.
	maxsites\tINT\tMaximum number of jackknife samples used for estimating std. err. (which sites are removed is random). Defaults to total number of available sites.
	nthr\tINT\tNumber of threads to use for estimating std. err. Default: 1
	oneCns\t0|1\t0: Only do Two-cns method. 1: Do both One-cns and Two-cns. Default: 0
	Arguments should be passed as: arg=value. Make sure that there is no spaces between arg, =, and value.
	Example call: Rscript ContaEst.R counts=example_counts freqs=HapMap.gz maxsites=500 nthr=20 outfile=outname oneCns=1
")

#fill hash with NAs
ArgHash<<-new.env()
for(i in PossibleArgNames){
	ArgHash[[i]]<-NA
}

#default values for bsize, oneCns and nthreads
ArgHash[["nthr"]]<-1
ArgHash[["oneCns"]]<-0
ArgHash[["maxsites"]]<-1e100 #make the default something huge in order to use all unless something else is specified

#remove spaces from argument names and vals
argnames<-gsub(" ", "", unlist(lapply(strsplit(args, "="), "[[", 1)))
#message(argnames)
argvals<-gsub(" ", "", unlist(lapply(strsplit(args, "="), "[[", 2)))
#message(argvals)

#check that all argnames are recognized
if(sum(!(argnames%in%PossibleArgNames))>0){
	message(USAGE)
	stop(paste(collapse=" ", c("Unrecognized arguments", argnames[!(argnames%in%PossibleArgNames)])))
}

#assign each arg val to an arg name in the arg hash
for(i in seq_along(argnames)){
	ArgHash[[argnames[i]]]<-argvals[i]
}

#retrieve arg vals from arghash, so they are used in an indep variable
countsname<-ArgHash[["counts"]]
contafrqsname<-ArgHash[["freqs"]]
outfilename<-ArgHash[["outfile"]]
nthr<-as.integer(ArgHash[["nthr"]])
doOneCns<-as.integer(ArgHash[["oneCns"]])
if(ArgHash[["maxsites"]]==1e100){
	maxnsites<-ArgHash[["maxsites"]]
}else{
	maxnsites<-as.integer(ArgHash[["maxsites"]])
}

#die if doOneCns is not 0 or 1
if(doOneCns!=0 & doOneCns!=1){
	message(USAGE)
	stop("oneCns can only be 0 or 1")	
}


#die if required arguments are missing
if(is.na(countsname) | is.na(contafrqsname) | is.na(outfilename)){
	message(USAGE)
	stop("counts, freqs and outfile are required")
}

#################

#args for testing
#countsname<-"/willerslev/scratch/bkl835/Contamination/ReWrite_250518/Tests/Yoruba_French_0.5_0.15_4C86u9.bam_CEU.icnts.gz_counts"
#contafrqsname<-"/willerslev/users-shared/science-snm-willerslev-bkl835/jmoreno/#ContaminationPanels/HapMapCEU.gz"
#outfilename<-"testout"
#maxnsites<-1000
#nthr<-30

require(doParallel)
registerDoParallel(cores=nthr)

bases<<-c("A", "C", "G", "T")

message("Reading ref panel... ")
freqs1<-read.table(gzfile(contafrqsname), as.is=T)

freqs1<-freqs1[freqs1[,4]=="+", ]

##############
#Run two-cns

freqs<-freqs1
colnames(freqs)<-NULL
rownames(freqs)<-NULL

#

message("Reading counts... ")
counts<-read.table(pipe(paste0("grep 'counts' ", countsname)), as.is=T)

# change to absolute positions:
counts[,2]<-counts[,2]+1

message("Calling cns... ")
# A C G T
GetCNS<-function(x){
	maxf<-as.numeric(which.max(x))
	if(sum(x==as.numeric(x[maxf]))==1){
		return(maxf)
	}else{
		rand<-sample(seq(1:sum(x==as.numeric(x[maxf]))), 1)
		return(which(x==as.numeric(x[maxf]))[rand])
		#return("N")
	}
}

cns<-apply(counts[,4:7], 1, GetCNS)

cns<-ifelse(cns==1, "A", ifelse(cns==2, "C", ifelse(cns==3, "G", "T")))

# Counts with all sites
newcounts<-cbind.data.frame(counts, cns, stringsAsFactors=F)

counts<-newcounts[,c(2, 4:8)]

message("Some extra filtering... ")
# change to absolute positions
freqs[,1]<-freqs[,1]+1

#create data.frame with sites
sitesC<-counts[counts[,1] %in% freqs[,1],]
sitesF<-freqs[freqs[,1] %in% sitesC[,1],c(2,5,3)]

sites<-data.frame("pos"=sitesC[,1], "A"=sitesC[,2], "C"=sitesC[,3], "G"=sitesC[,4], "T"=sitesC[,5], "cns"=sitesC[,6], "A1_1"=sitesF[,1], "A2_1"=sitesF[,2], "frq1"=sitesF[,3],  stringsAsFactors=F)

#columns in sites:
#1. pos
#2. a counts
#3. c counts
#4. g counts
#5. t counts
#6. consensus
#7. ref allele
#8. alt allele in contaminant
#9. freq of ref allele in contaminant population

# remove the sites that are within 5 bases from each other
DistS<-abs(c(sites$pos[1:(length(sites$pos)-1)]-sites$pos[2:length(sites$pos)], 1000))
rem<-unique(c(which(DistS<=5), which(DistS<=5)+1))

if(length(rem)>0){
	sites<-sites[-rem,]
}

#remove triallelic sites
sites<-sites[(sites$cns==sites$A1_1 | sites$cns==sites$A2_1), ]

message("Computing overall error rate... ")
###
##This is for when the c pgm does the filtering
#get counts for adjacent positions
AdjPos<-unique(c(sites$pos-1, sites$pos-2, sites$pos-3, sites$pos-4, sites$pos-5, sites$pos+1, sites$pos+2, sites$pos+3, sites$pos+4, sites$pos+5))
AdjCnts<-counts[counts[,1] %in% AdjPos, ]
###

####

#get overall error rate

cntMisMatches<-function(x){
	x<-x[-1]
	if(x[5]=="A"){
		mm<-sum(as.numeric(x[-5][-1]))
	}
	if(x[5]=="C"){
		mm<-sum(as.numeric(x[-5][-2]))
	}
	if(x[5]=="G"){
		mm<-sum(as.numeric(x[-5][-3]))
	}
	if(x[5]=="T"){
		mm<-sum(as.numeric(x[-5][-4]))
	}
	return(mm)
}

###
##This is for when the c pgm does the filtering
MMcnts<-unname(apply(AdjCnts, 1, cntMisMatches))
AdjLens<-unname(rowSums(AdjCnts[,2:5]))
###

MMrate<<-sum(MMcnts)/sum(AdjLens)
sites<<-sites

#two-consensus likelihood

ParseRow<-function(x){
	cnsbase<-x[6]
	cnsct<-as.numeric(x[which(bases==cnsbase)+1])
	altbase<-x[c(7,8)][x[c(7,8)]!=cnsbase]
	altct<-sum(as.numeric(x[which(bases!=cnsbase)+1]))
	if(cnsbase==x[7]){
		altfrq<-(1-as.numeric(x[9]))
	}else if(cnsbase==x[8]){
		altfrq<-as.numeric(x[9])
	}
	return(c(cnsct, altct, altfrq))
}

cnstable<<-t(apply(sites, 1, ParseRow))

###
#two-consensus likelihood

message("Running two-cns... ")

TwoConsmLL<-function(conta){
	fis<-(conta*cnstable[,3]*(1-((4/3)* MMrate)))+MMrate
	fisinv<-(conta*(1-cnstable[,3])*(1-((4/3)* MMrate)))+MMrate
	probs<-(choose((cnstable[,1]+cnstable[,2]), cnstable[,2])*(fis^cnstable[,2])*((1-fis)^cnstable[,1]))/2
	probsinv<-(choose((cnstable[,1]+cnstable[,2]), cnstable[,1])*(fisinv^cnstable[,1])*((1-fisinv)^cnstable[,2]))/2
	probs<-log(probs+probsinv)
	return(-sum(probs))
}

est2<-optim(0.25, TwoConsmLL, lower=c(0.000000000001), upper=c(1-0.000000000001), method="Brent")$par

message("Jackknife for two-cns... ")
#get s.e. and 95 CI

#decide which sites to keep
nsites<-nrow(cnstable)
if(nsites > maxnsites){
	nsites<-maxnsites
	sitestorm<-sample(1:nrow(cnstable), maxnsites, replace=F)
}else{
	sitestorm<-1:nrow(cnstable)
}

#jk estimate array

jke<-unlist(foreach(i=seq_along(sitestorm)) %dopar%{
	currcnstable<-cnstable[-sitestorm[i], ]
	TwoConsmLL<-function(conta){
		fis<-(conta*currcnstable[,3]*(1-((4/3)* MMrate)))+MMrate
		fisinv<-(conta*(1-currcnstable[,3])*(1-((4/3)* MMrate)))+MMrate
		probs<-(choose((currcnstable[,1]+currcnstable[,2]), currcnstable[,2])*(fis^currcnstable[,2])*((1-fis)^currcnstable[,1]))/2
		probsinv<-(choose((currcnstable[,1]+currcnstable[,2]), currcnstable[,1])*(fisinv^currcnstable[,1])*((1-fisinv)^currcnstable[,2]))/2
		probs<-log(probs+probsinv)
		return(-sum(probs))
	}
	optim(0.25, TwoConsmLL, lower=c(0.000000000001), upper=c(1-0.000000000001), method="Brent")$par
})

#get sigma and upper & lower bounds for 95%CI

s<-sqrt(((nsites-1)/nsites)*(sum((jke-est2)^2)))

upper2<-est2+qnorm(1-(.05)/2)*s
lower2<-est2-qnorm(1-(.05)/2)*s

##################
#one-cns method, only if the user requires it

if(doOneCns==1){
#one-consensus likelihood

message("Running one-cns... ")

ConsmLL<-function(conta){
	fis<-(conta*cnstable[,3]*(1-((4/3)*MMrate)))+MMrate
	probs<-log(choose((cnstable[,1]+cnstable[,2]), cnstable[,2])*(fis^cnstable[,2])*((1-fis)^cnstable[,1]))
	return(-sum(probs))
}

est1<-optim(0.25, ConsmLL, lower=c(0.000000000001), upper=c(1-0.000000000001), method="Brent")$par

message("Jackknife for one-cns... ")
#get s.e. and 95 CI

#jk estimate array

jke<-unlist(foreach(i=seq_along(sitestorm)) %dopar%{
	currcnstable<-cnstable[-sitestorm[i], ]
	ConsmLL<-function(conta){
		fis<-(conta*currcnstable[,3]*(1-((4/3)* MMrate)))+MMrate
		probs<-log(choose((currcnstable[,1]+currcnstable[,2]), currcnstable[,2])*(fis^currcnstable[,2])*((1-fis)^currcnstable[,1]))
		return(-sum(probs))
	}
	optim(0.25, ConsmLL, lower=c(0.000000000001), upper=c(1-0.000000000001), method="Brent")$par
})


#get sigma and upper & lower bounds for 95%CI

s<-sqrt(((nsites-1)/nsites)*(sum((jke-est1)^2)))

upper1<-est1+qnorm(1-(.05)/2)*s
lower1<-est1-qnorm(1-(.05)/2)*s
}
##################

if(doOneCns==0){
	writeLines(paste(sep="\t", "Two-cns", est2, lower2, upper2, MMrate, nrow(cnstable)), con= outfilename)
}else{
	writeLines(c(paste(sep="\t", "One-cns", est1, lower1, upper1, MMrate, nrow(cnstable)), paste(sep="\t", "Two-cns", est2, lower2, upper2, MMrate, nrow(cnstable))), con= outfilename)
}









