/*
  This is a class that dumps plink file output

 */
#include <assert.h>

#include "analysisFunction.h"
#include "shared.h"
#include <htslib/kstring.h>
#include "abcCallGenotypes.h"
#include "abcWritePlink.h"
void abcWritePlink::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doPlink\t%d\n",doPlink);
  fprintf(argFile,"\t1: binary fam/bim/bed format (still beta, not really working)\n");
  fprintf(argFile,"\t2: tfam/tped format\n");
  fprintf(argFile,"\n\tNB This is a wrapper around -doGeno see more information for that option\n");
}

void abcWritePlink::run(funkyPars *pars){
  if(doPlink==0)
    return ;
   
}

void abcWritePlink::clean(funkyPars *pars){
  if(doPlink==0)
    return;
    

}

unsigned char recode(int i){
  if(i==2)
    return '\x03';
  else if(i==1)
    return '\x02';
  else if(i==0)
    return '\x00';
  else if(i==-1)
    return '\x01';
  assert(0==1);  
  return 'x';//to remove compile warnings;
}

void abcWritePlink::print(funkyPars *pars){
  if(doPlink==0)
    return;
  genoCalls *geno =(genoCalls *) pars->extras[10];
  char *chr = header->target_name[pars->refId];
  if(doPlink==2){
    kstring_t kstr;
    kstr.s=NULL;kstr.l=kstr.m=0;
    
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      ksprintf(&kstr,"%s %s_%d 0 %d",chr,chr,pars->posi[s]+1,pars->posi[s]+1);
      
      for(int i=0;i<pars->nInd;i++){
      	if(geno->dat[s][i]==0)
	  ksprintf(&kstr,"\t%c %c",intToRef[pars->major[s]],intToRef[pars->major[s]]);
	else if(geno->dat[s][i]==1)
	  ksprintf(&kstr,"\t%c %c",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
	else if(geno->dat[s][i]==2)
	  ksprintf(&kstr,"\t%c %c",intToRef[pars->minor[s]],intToRef[pars->minor[s]]);
	else if(geno->dat[s][i]==-1)
	  ksprintf(&kstr,"\t0 0");
      }
      kputc('\n',&kstr);
    }
    fwrite(kstr.s,1,kstr.l,fp2);
    free(kstr.s);
  }else if(doPlink==1){

    kstring_t kstr;
    kstr.s=NULL;kstr.l=kstr.m=0;


    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
     
      ksprintf(&kstr,"%s\t%s_%d\t0\t%d",chr,chr,pars->posi[s]+1,pars->posi[s]+1); 
      
      int minAl=0;
      int majAl=0;
      for(int i=0;i<pars->nInd;i++)
	if(geno->dat[s][i]!=-1){
	  minAl+=geno->dat[s][i];
	  majAl+=2-geno->dat[s][i];
	}
      fprintf(stderr,"posi:%d: maj:%d min:%d\n",pars->posi[s]+1,minAl,majAl);
      int swapped =0;
      if(majAl!=0&&minAl!=0){
	if(majAl>minAl)
	  ksprintf(&kstr,"\t%c\t%c\n",intToRef[pars->minor[s]],intToRef[pars->major[s]]);
	else if(majAl<minAl){
	  ksprintf(&kstr,"\t%c\t%c\n",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
	  swapped =1;
	}else if(majAl==minAl){
	  int mmin=std::min(pars->major[s],pars->minor[s]);
	  if(mmin==pars->major[s])
	    ksprintf(&kstr,"\t%c\t%c\n",intToRef[pars->minor[s]],intToRef[pars->major[s]]);
	  else{
	    ksprintf(&kstr,"\t%c\t%c\n",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
	    swapped =1;
	  }
	}
      }else{
	ksprintf(&kstr,"\t%c\t%c\n",intToRef[pars->minor[s]],intToRef[pars->major[s]]);

      }
      
      
     

      unsigned char part;
      int idx =0;

      
      
      for(int i=0;i<pars->nInd;i++){
	unsigned char obs = recode(geno->dat[s][i]);
	//	fprintf(stderr,"%d:%d\n",geno->dat[s][i],(int)obs);
	if(idx==4){
	  fwrite(&part,1,1,fp2);
	  idx=0;
	}
	if(idx==0)
	  part = obs;
	else
	  part |= obs<<(2*idx);
	idx++;
      }
      if(idx)
	fwrite(&part,1,1,fp2);
      
    }
    fwrite(kstr.s,1,kstr.l,fp1);
    free(kstr.s);



  }
  
}

void abcWritePlink::getOptions(argStruct *arguments){

  doPlink=angsd::getArg("-doPlink",doPlink,arguments);
  if(doPlink==0)
    return;
  
  int doGeno = 0;
  doGeno = angsd::getArg("-doGeno",doGeno,arguments);

  if(doPlink==1){
    fprintf(stderr,"Binary plink is not suported yet. use -doPlink 2 \n");
    exit(0);

  }

  if(doPlink!=0 && doGeno==0){
    fprintf(stderr,"Must supply -doGeno  to write plink files (consider supplying negative value for suprresing .geno.gz output)\n");
    exit(0);
  }
  assert(doPlink>=0 && doPlink<=2);
  if(doPlink==0)
    return;
  

}

void writeTfam(FILE *fp,int nInd){
  for(int i=1;i<=nInd;i++)
    fprintf(fp,"%d 1 0 0 0 -9\n",i);
  
}



abcWritePlink::abcWritePlink(const char *outfiles,argStruct *arguments,int inputtype){
  fp1=fp2=NULL;
  doPlink =0;
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doPlink")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);


  if(doPlink==0){
    shouldRun[index] =0;
    return;
}
  printArg(arguments->argumentFile);
  if(doPlink == 1){
    const char* postfix1=".fam";
    fp1 = aio::openFile(outfiles,postfix1);
    const char* postfix2=".bed";
    fp2 = aio::openFile(outfiles,postfix2);
    //write magic nubmers
    fputc('\x6C',fp2);
    fputc('\x1B',fp2);
    fputc('\x01',fp2);
  }else if (doPlink==2){
    const char* postfix1=".tfam";
    fp1 = aio::openFile(outfiles,postfix1);
    const char* postfix2=".tped";
    fp2 = aio::openFile(outfiles,postfix2);
  }
  writeTfam(fp1,arguments->nInd);
  if(fp1) fclose(fp1);fp1=NULL;
  if(doPlink==1){
    const char* postfix2=".bim";
    fp1 = aio::openFile(outfiles,postfix2);
  }
  
  
}


abcWritePlink::~abcWritePlink(){
  if(fp1)     fclose(fp1);
  if(fp2)     fclose(fp2);
}


