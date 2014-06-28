/*
  This is a modified version of the code recieved by EUNJUNG HAN from jnovembre lab.
This is a banded version of the dynamic algorithm used for calculating the sample allele frequncy, or multisample GLs or site likelihoods.


 */

#include <cstdio>
#include <cmath>
#include <assert.h>
#include <zlib.h>
#include <algorithm>
#include "kstring.h"
#include "abcFreq.h"
#include "shared.h"
#include "analysisFunction.h"

#include "abc.h"

#include "abcSaf2.h"


#define MINLIKE -1000.0 //this is for setting genotypelikelhoods to missing (EXPLAINED BELOW)


double lbico(double n, double k);


void abcSaf2::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSaf2\t\t%d\n",doSaf2);
  fprintf(argFile,"\t1: perform multisample GL estimation EUNJUNG VERSION(cite that, more info )\n");
  fprintf(argFile,"\t-underFlowProtect\t%d\n",underFlowProtect);
  fprintf(argFile,"\t-fold\t\t\t%d (deprecated)\n",fold);
  fprintf(argFile,"\t-anc\t\t\t%s (ancestral fasta)\n",anc);
  fprintf(argFile,"\t-noTrans\t\t%d (remove transitions)\n",noTrans);
  fprintf(argFile,"\t-mode\t\t\t%d (Dynamic programming algorithm to compute site likelihood vectors)\n\t0: original dynamic programming (default)\n\t1: rescaled dynamic programming\n\t2: adaptive K-restricted dynamic programming\n",mode); //[EJ]
}


void abcSaf2::getOptions(argStruct *arguments){
  doSaf2=angsd::getArg("-doSaf2",doSaf2,arguments);

  noTrans = angsd::getArg("-noTrans",noTrans,arguments);
  mode = angsd::getArg("-mode",mode,arguments); //[EJ]
    
  int GL = 0;
  GL = angsd::getArg("-GL",GL,arguments);
  
  underFlowProtect=angsd::getArg("-underFlowProtect",underFlowProtect,arguments);
  fold=angsd::getArg("-fold",fold,arguments);
  char *sim1 = NULL;
  isSim=angsd::getArg("-isSim",isSim,arguments);

  //  fprintf(stderr,"sim1=%p %s\n",sim1,sim1);

  if(doSaf2==0)
    return;
  anc = angsd::getArg("-anc",anc,arguments);
  if(doSaf2 && (anc==NULL&&isSim==0)){
    if(doSaf2!=3){
      fprintf(stderr,"Must supply -anc for polarizing the spectrum\n");
      exit(0);
    }
  }

  if(GL==0 &&arguments->inputtype==7){//DRAGON this is not enough, but this is most likely what everyone done...
    fprintf(stderr,"Must supply genotype likelihoods (-GL [INT])\n");
    printArg(arguments->argumentFile);
    exit(0);
  }

}

abcSaf2::abcSaf2(const char *outfiles,argStruct *arguments,int inputtype){
  const char *SAF = ".saf";
  const char *SFSPOS =".saf.pos.gz";

  //default
  underFlowProtect = 0;
  fold =0;
  isSim =0;
  //from command line
  anc=NULL;
  
  noTrans = 0;
  
  mode=0; //0,1,2 [EJ]
  doSaf2=0;
  
  outfileSFS = NULL;
  outfileSFSPOS = Z_NULL;
  
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSaf2")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }



  getOptions(arguments);

  printArg(arguments->argumentFile);  
  if(doSaf2==0){
    shouldRun[index] =0;
    return;
  }
  newDim = 2*arguments->nInd+1;
  if(fold)
    newDim = arguments->nInd+1;

  if(doSaf2==1){
    outfileSFS =  aio::openFile(outfiles,SAF);
    outfileSFSPOS =  aio::openFileGz(outfiles,SFSPOS,GZOPT);
  }
   
}


abcSaf2::~abcSaf2(){

  if(outfileSFS) fclose(outfileSFS);;
  if(outfileSFSPOS) gzclose(outfileSFSPOS);
}

void abcSaf2::run(funkyPars  *p){
  if(p->numSites==0||(doSaf2==0 ))
    return;
  if(doSaf2!=1){
    fprintf(stderr,"\t-> Error -doSaf22 is only implemented for simple SFS estimation\n");
    exit(0);
  }
  
  realRes2 *r = new realRes2;
  r->oklist=new char[p->numSites];
  memset(r->oklist,0,p->numSites);
  r->pLikes=new double*[p->numSites];
  p->extras[index] = r;
  
  algoJoint(p->likes,p->anc,p->numSites,p->nInd,underFlowProtect,fold,p->keepSites,r,noTrans);
  
}

void abcSaf2::clean(funkyPars *p){
  if(p->numSites==0||(doSaf2==0 ))
    return;
  if(doSaf2==3)
    return;
  realRes2 *r=(realRes2 *) p->extras[index];
  
  //  realRes2 *r=(realRes2 *) p->extras[index];
  int id=0;
  for(int i=0;i<p->numSites;i++)
    if(r->oklist[i]==1)
      delete [] r->pLikes[id++];
  delete [] r->pLikes;
  delete [] r->oklist;
  delete r;

}

void abcSaf2::print(funkyPars *p){
  if(p->numSites==0||(doSaf2==0))
    return;
  //  fprintf(stderr,"newDim:%d doSaf2:%d\n",newDim,doSaf2);
  if(doSaf2==1){
    realRes2 *r=(realRes2 *) p->extras[index];
    int id=0;
    void printFull(funkyPars *p,int index,FILE *outfileSFS,gzFile outfileSFSPOS,char *chr,int newDim);
    printFull(p,index,outfileSFS,outfileSFSPOS,header->name[p->refId],newDim);
  }   
}


//---------------------------------------------------------------------------------------------------------------------------//
// Added my EJ
inline double getLikelihoods(double *L012, double *like, int i, int AA_offset, int Aa_offset, int aa_offset,int underFlowProtect) __attribute__((always_inline));
 double getLikelihoods(double *L012, double *like, int i, int AA_offset, int Aa_offset, int aa_offset,int underFlowProtect){
    double GAA,GAa,Gaa,mymax;
    
    // genotye likelihhods (het prob is multiplied by 2) in log space
    GAA = like[i*10+AA_offset];
    GAa = like[i*10+Aa_offset]+log(2.0);
    Gaa = like[i*10+aa_offset];
    
    // normalize
    if (Gaa > GAa && Gaa > GAA) mymax = Gaa;
    else if (GAa > GAA)         mymax = GAa;
    else                        mymax = GAA;
    
    if(mymax<MINLIKE){
        Gaa = 0;
        GAa = 0;
        GAA = 0;
        //totmax = totmax + mymax;
    }else{
        Gaa=Gaa-mymax;
        GAa=GAa-mymax;
        GAA=GAA-mymax;
        //totmax = totmax + mymax;
    }
    
    if(underFlowProtect==0){//no protection
      L012[0]=exp(Gaa);      L012[1]=exp(GAa);      L012[2]=exp(GAA);
    }
    else{
      L012[0]=Gaa;      L012[1]=GAa;      L012[2]=GAA;
    }
    //fprintf(stderr,"aa=%d,Aa=%d,AA=%d: L0=%f\tL1=%f\tL2=%f\n",aa_offset,Aa_offset,AA_offset,L0,L1,L2);
    return mymax;
}


double findMax(double *array, int start, int end){
  assert(end-start>0);
  double mmax=array[start];
  for(int i=start+1;i<end;i++)
    if(array[i]>mmax)
      mmax=array[i];
  return mmax;
}

void dynamic(double *like, double *sumMinors, int nInd,int AA_offset, int Aa_offset, int aa_offset, int underFlowProtect){
    int nChr=2*nInd;
    int nBins=nChr+1;
    double mymax,totmax;
    double L012[3];
    double hj[nBins];
    for(int i=0;i<nBins;i++) hj[i]=.0;
    totmax=0.0;
    
    // Update hj for j=0
    mymax=getLikelihoods(hj,like,0,AA_offset,Aa_offset,aa_offset,underFlowProtect);
    totmax+=mymax;
    
    // Update hj for j>0
    for(int j=1 ; j<nInd ;j++) {//for each individual...
        mymax=getLikelihoods(L012,like,j,AA_offset,Aa_offset,aa_offset,underFlowProtect);
        totmax+=mymax;
        for(int k=2*(j+1); k>1;k--){
            hj[k] = L012[2]*hj[k-2]+L012[1]*hj[k-1]+L012[0]*hj[k];
	    assert(!std::isnan(hj[k]));
        }
        hj[1] = L012[0]*hj[1]+L012[1]*hj[0];
        hj[0] = L012[0]*hj[0];
    } //end of individulas; if underflow protect then hj is in logspace
    
    // Normalize hj - log(hj[i])=log(Zn[i]/(2n,i)) + Add it to sumMinors (in normal space)
    for(int i=0;i<nBins;i++) sumMinors[i] += exp(log(hj[i])-lbico(nChr,i)+totmax);
}

void dynamicAdaptiveK(double *like, double *sumMinors, int nInd,int AA_offset, int Aa_offset, int aa_offset, int underFlowProtect){

    double epsilon=0.000001; //10^(-6)
    int nchrom=2*nInd;
    int nBins=2*nInd+1;
    int from, fromOld,From, to, GT, nChr;
    double mymax,totmax,checkVal,l1;
    double L012[3];
    // Initialize hj
    double hj[nBins];
    for(int i=0;i<nBins;i++) hj[i]=.0;
    totmax=0.0;
    
    // Update hj for j=0
    mymax=getLikelihoods(hj,like,0,AA_offset,Aa_offset,aa_offset,underFlowProtect);
    totmax+=mymax;

    from=0; to=2;
    
    // Update hj for i>0
    for(int j=1 ; j<nInd ;j++) {
        mymax=getLikelihoods(L012,like,j,AA_offset,Aa_offset,aa_offset,underFlowProtect);
        totmax+=mymax;
        nChr=2*(j+1);
     	GT=0; 
	l1=L012[1]/2.0;
        if(l1>L012[0] && l1>L012[2]) GT=1;
        else if(L012[2]>L012[0] && L012[2]>l1) GT=2;
  
        
        // range
        fromOld=from;
        from+=GT; to+=GT;
        
        // Check point (left boundary)
        if(from>0){
            if(from>1) checkVal = L012[2]*hj[from-2]+L012[1]*hj[from-1]+L012[0]*hj[from];
            else       checkVal = L012[0]*hj[1]+L012[1]*hj[0];
            
            while(checkVal>epsilon) {
                from--;
                if(from==0) break;
                if(from>1) checkVal = L012[1]*hj[from-2]+L012[1]*hj[from-1]+L012[0]*hj[from];
                else       checkVal = L012[0]*hj[1]+L012[1]*hj[0];
            }
        }
        
        // Check point (right boundary)
        if(to<nChr){
            checkVal = L012[2]*hj[to-2]+L012[1]*hj[to-1]+L012[0]*hj[to];
            while(checkVal>epsilon) {
                to++;
                if(to==nChr) break;
                checkVal = L012[2]*hj[to-2]+L012[1]*hj[to-1]+L012[0]*hj[to];            
	    }
        }
	
        //update
        From=fmax(2,from);
        for(int k=to; k>=From;k--){
            hj[k] = L012[2]*hj[k-2]+L012[1]*hj[k-1]+L012[0]*hj[k];
            assert(!std::isnan(hj[k]));
        }
        if(from <= 1) {
            hj[1]= L012[0]*hj[1]+L012[1]*hj[0];
            if(from==0) hj[0]=L012[0]*hj[0];
        }
        while(fromOld<from){
            hj[fromOld]=0;
            fromOld++;
        }
    }
    
    // Add it to sumMinors (in normal space)
    for(int k=from; k<=to; k++) sumMinors[k] += exp(log(hj[k])-lbico(nchrom,k)+totmax);
}



void rescaledDynamic(double *like, double *sumMinors, int nInd,int AA_offset, int Aa_offset, int aa_offset, int underFlowProtect){
    int nChrom=2*nInd;
    int nBins=nChrom+1;
    int nChr;
    double mymax,totmax, maxhj,denom,rescalingFactor;

    double L012[3];
    // Initialize hj
    double hj[nBins];
    for(int i=0;i<nBins;i++) hj[i]=.0;
    totmax=0.0;

    // Update hj for j=0
    mymax=getLikelihoods(hj,like,0,AA_offset,Aa_offset,aa_offset,underFlowProtect);
    hj[1] /=2.0;
    totmax+=mymax;
    nChr=2;
    
    // Update hj for i>0
    rescalingFactor=0;
    for(int j=1 ; j<nInd ;j++) {
        //rescaling
        maxhj=findMax(hj,0,nChr);
        rescalingFactor+=log(maxhj);

	for(int i=0;i<=nChr;i++)
	  hj[i] /= maxhj;

        // update
        mymax=getLikelihoods(L012,like,j,AA_offset,Aa_offset,aa_offset,underFlowProtect);

        totmax+=mymax;
        nChr=2*(j+1); denom=1.0/(nChr*(nChr-1));
        for(int k=nChr; k>1;k--){
	  hj[k] = ((nChr-k)*(nChr-k-1)*L012[0]*hj[k] + k*(nChr-k)*L012[1]*hj[k-1]+k*(k-1)*L012[2]*hj[k-2])*denom;
	  assert(!std::isnan(hj[k]));
        }
        hj[1] = ((nChr-2)*L012[0]*hj[1]+L012[1]*hj[0])/nChr;
        hj[0] = L012[0]*hj[0];
    }

     // Rescale hj + Add it to sumMinors (in normal space)
    for(int k=0; k<nBins; k++)
        sumMinors[k] += exp(log(hj[k])+totmax+rescalingFactor);
}


void adaptiveK(double *like, double *sumMinors, int nInd,int AA_offset, int Aa_offset, int aa_offset, int underFlowProtect){
    //fprintf(stderr,"Adaptive K algorithm\n");
    double epsilon=1e-6; //10^(-6)
    int nBins=2*nInd+1;
    int from, fromOld,From, to, GT, nChr;
    double mymax,totmax, maxhj,checkVal, denom,l1, rescalingFactor;
    double L012[3];
    // Initialize hj
    double hj[nBins];    for(int i=0;i<nBins;i++) hj[i]=.0;
    totmax=0.0;
    
    // Update hj for j=0
    mymax=getLikelihoods(hj,like,0,AA_offset,Aa_offset,aa_offset,underFlowProtect);
    hj[1] /= 2.0;
    totmax+=mymax;
    nChr=2;
    from=0; to=2;
    
    // Update hj for i>0
    rescalingFactor=0;
    for(int j=1 ; j<nInd ;j++) {
        //rescaling
        maxhj=findMax(hj,from,to);
        rescalingFactor+=log(maxhj);
	for(int i=0;i<=nChr;i++)
	  hj[i] /= maxhj;

        //update
        mymax=getLikelihoods(L012,like,j,AA_offset,Aa_offset,aa_offset,underFlowProtect);
        totmax+=mymax;
        nChr=2*(j+1); denom=1.0/(nChr*(nChr-1));
        
	GT=0; 
	l1=L012[1]/2.0;
        if(l1>L012[0] && l1>L012[2]) GT=1;
        else if(L012[2]>L012[0] && L012[2]>l1) GT=2;
        
        // range
        fromOld=from;
        from+=GT; to+=GT;
        
        // Check point (left boundary)
        if(from>0){
            if(from>1)
                checkVal = ((nChr-from)*(nChr-from-1)*L012[0]*hj[from]+from*(nChr-from)*L012[1]*hj[from-1]+from*(from-1)*L012[2]*hj[from-2])*denom;
            else
                checkVal = ((nChr-2)*L012[0]*hj[1]+L012[1]*hj[0])/nChr;
                
            while(checkVal>epsilon) {
                    from--;
                    if(from==0) break;
                    if(from>1)
                        checkVal = ((nChr-from)*(nChr-from-1)*L012[0]*hj[from]+from*(nChr-from)*L012[1]*hj[from-1]+from*(from-1)*L012[2]*hj[from-2])*denom;
                    else
                        checkVal = ((nChr-2)*L012[0]*hj[1]+L012[1]*hj[0])/nChr;
            }
        }
        
        // Check point (right boundary)
        if(to<nChr){
            checkVal = ((nChr-to)*(nChr-to-1)*L012[0]*hj[to]+to*(nChr-to)*L012[1]*hj[to-1]+to*(to-1)*L012[2]*hj[to-2])*denom;
            while(checkVal>epsilon) {
                to++;
                if(to==nChr) break;
                checkVal = ((nChr-to)*(nChr-to-1)*L012[0]*hj[to]+to*(nChr-to)*L012[1]*hj[to-1]+to*(to-1)*L012[2]*hj[to-2])*denom;
            }
        }
        //fprintf(stderr,"%d: %d - %d\n",j+1,from,to);
        
        //update
        From=fmax(2,from);
        for(int k=to; k>=From;k--){
	  hj[k] = denom* ( (nChr-k)*(nChr-k-1)*L012[0]*hj[k]+k*(nChr-k)*L012[1]*hj[k-1]+k*(k-1)*L012[2]*hj[k-2]);
          assert(!std::isnan(hj[k]));
        }
        if(from <= 1) {
            hj[1]=((nChr-2)*L012[0]*hj[1]+L012[1]*hj[0])/nChr;
            if(from==0) hj[0]=L012[0]*hj[0];
        }
        while(fromOld<from){
            hj[fromOld]=0;
            fromOld++;
        }
    }
    
    // Rescale hj + Add it to sumMinors (in normal space)
    for(int k=from; k<=to; k++)
        sumMinors[k] += exp(log(hj[k])+totmax+rescalingFactor);
}

void abcSaf2::algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes2 *r,int noTrans) {
    // fprintf(stderr,"[%s]\n",__FUNCTION__);
    int myCounter =0;
    if(anc==NULL||liks==NULL){
        fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n",__FUNCTION__,liks,anc);
        exit(0);
    }
    double sumMinors[2*numInds+1];  //the sum of the 3 different minors
    
    for(int it=0; it<nsites; it++) {//loop over sites
        double *like = liks[it]; //[EJ]
        int major_offset = anc[it];
        if(major_offset==4||(keepSites[it]==0)){//skip of no ancestral information
            //      r->oklist is zero no need to update
            continue;
        }
        //set the resultarray to zeros
        for(int sm=0 ; sm<(2*numInds+1) ; sm++ )
            sumMinors[sm] = 0;
        
        //loop through the 3 different minors
        for(int minor_offset=0;minor_offset<4;minor_offset++) {
            if(minor_offset == major_offset)
                continue;
            if(noTrans){
                if((major_offset==2&&minor_offset==0)||(major_offset==0&&minor_offset==2))
                    continue;
                if((major_offset==1&&minor_offset==3)||(major_offset==3&&minor_offset==1))
                    continue;
            }
           
            //hook for only calculating one minor
            int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
            int AA_offset = angsd::majorminor[minor_offset][minor_offset];//0-9
            int aa_offset = angsd::majorminor[major_offset][major_offset];//0-9
            //      fprintf(stderr,"%d:%d\t%d\t%d\n",major_offset,Aa_offset,AA_offset,aa_offset);
            
            //--------------part two ------------------------//
            if(mode==0)      dynamic(like,sumMinors,numInds,AA_offset,Aa_offset,aa_offset,underFlowProtect);          //[EJ]
	    else if(mode==1) rescaledDynamic(like,sumMinors,numInds,AA_offset,Aa_offset,aa_offset,underFlowProtect);  //[EJ]

	    else if(mode==2) adaptiveK(like,sumMinors,numInds,AA_offset,Aa_offset,aa_offset,underFlowProtect);        //[EJ]
            else if(mode==3) dynamicAdaptiveK(like,sumMinors,numInds,AA_offset,Aa_offset,aa_offset,underFlowProtect);        //[EJ]
            //------------- end of part two -----------------------------------------------//
#if 0
            for(int ii=0;ii<10*numInds;ii++) fprintf(stdout,"%f\t",liks[it][ii]);
#endif
        }

        //sumMinors is in normal space, not log
        /*
         we do 3 things.
         1. log scaling everyting
         2. rescaling to the most likely in order to avoid underflows in the optimization
         3. we might do a fold also.
         
         */
        
        if(fold) {
            int newDim = numInds+1;
            for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
                sumMinors[i] = log(sumMinors[i] + sumMinors[2*numInds-i]);//THORFINN NEW
            sumMinors[newDim-1] = log(sumMinors[newDim-1])+log(2.0);
            angsd::logrescale(sumMinors,newDim);
            if(std::isnan(sumMinors[0]))
                r->oklist[it] = 2;
            else{
                r->oklist[it] = 1;
                r->pLikes[myCounter] =new double[numInds+1];
                memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(numInds+1));
                myCounter++;
            }
        }else{
            for(int i=0;i<2*numInds+1;i++)
                sumMinors[i] = log(sumMinors[i]);
            angsd::logrescale(sumMinors,2*numInds+1);
            if(std::isnan(sumMinors[0]))
                r->oklist[it] = 2;
            else{
                r->oklist[it] = 1;
                r->pLikes[myCounter] =new double[2*numInds+1];
                memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(2*numInds+1));
                myCounter++;
            }
        
	    #if 0
            // [EJ] Print output to STD ERR
            fprintf(stdout,"%d\t",it+1);//[EJ]
            for(int i=0;i<(2*numInds+1);i++) fprintf(stdout,"%f\t",(sumMinors[i])); //[EJ]
            fprintf(stdout,"\n"); //[EJ]
	    exit(0);
	    #endif
	}
    }
    
}

