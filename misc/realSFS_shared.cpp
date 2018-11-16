#include "realSFS_shared.h"
int **posiG  = NULL;

//unthreaded
//this will populate the keep vector by
// 1) set the chooseChr and populate toKeep
// 2) find over lap between different positions
// this is run once for each chromsome

int set_intersect_pos(std::vector<persaf *> &saf,char *chooseChr,int start,int stop,char **curChr,filt *fl){
  //fprintf(stderr,"[%s] chooseChr:%s, start:%d stop:%d\n",__FUNCTION__,chooseChr,start,stop );

  if(0&&saf.size()==1&&chooseChr==NULL){//use entire genome, then don't do any strange filtering
    //fprintf(stderr,"herer\n");
    return 0 ;
  }
  /*
    What happens here? very good question

   */
  static int firstTime =1;
  static myMap::iterator it_outer=saf[0]->mm.begin();
 aGotoHereIsTheEasiest:
  if(chooseChr==NULL){
    if(it_outer==saf[0]->mm.end())//if we are done
      return -2;
    else if(firstTime==0){
      it_outer++;
      if(it_outer==saf[0]->mm.end()){
	//	fprintf(stderr,"done reading will exit \n");
	return -3;
      }
    }else
      firstTime =0;
    chooseChr = it_outer->first;
    if(curChr)
      *curChr=it_outer->first;
    fprintf(stderr,"\t-> Is in multi sfs, will now read data from chr:%s\n",chooseChr);
  }
  fprintf(stderr,"\t-> hello Im the master merge part of realSFS. and I'll now do a tripple bypass to find intersect \n");
  fprintf(stderr,"\t-> 1) Will set iter according to chooseChr and start and stop, and possibly using -sites\n");
  assert(chooseChr!=NULL);
  //check that chooseChr exists in all pops.
  for(int j=0;j<saf.size();j++){
    myMap::iterator it = saf[j]->mm.find(chooseChr);
    if(it==saf[j]->mm.end()){
      fprintf(stderr,"\t-> Chromosome: \'%s\' does not exists in population: %s will skip it\n",chooseChr,saf[j]->fname);
      //fprintf(stderr,"asdf: %s\n",it_outer->first);
      chooseChr=NULL;
      goto aGotoHereIsTheEasiest;
    }

  }


  //hit will contain the depth across different populations
  keep<char> *hit =NULL;

  //  if(saf.size()>1)
  hit =keep_alloc<char>();//
  
  
  //this loop will populate a 'hit' array containing the effective (differnt pops) depth
  //if we only have one population, then just return after iter_init
  int killbreak =0;
  for(int i=0;i<saf.size();i++) {
    myMap::iterator it = iter_init(saf[i],chooseChr,start,stop);
    if(fl!=NULL&&i==0){
      filt_readSites(fl,chooseChr,saf[0]->ppos[it->second.nSites-1]);
      if(fl!=NULL&&fl->keeps==NULL){
	chooseChr=NULL;
	goto aGotoHereIsTheEasiest;
      }

    }
    
    //    fprintf(stderr,"ASDFASDF:%p\n",saf[i]->ppos);
    //    assert(it!=saf[i]->mm.end());  
    if(it==saf[i]->mm.end()){
      killbreak =1;
      break;
    }
    if(saf.size()==1 &&fl==NULL){
      keep_destroy(hit);
      return 0;
    }if(saf[i]->ppos[it->second.nSites-1] >= hit->m)
      realloc(hit,saf[i]->ppos[it->second.nSites-1]+1);
    assert(hit->m>0);
    for(size_t j=saf[i]->toKeep->first;j<=saf[i]->toKeep->last;j++){
      if(j>saf[i]->toKeep->first&&saf[i]->ppos[j] <= saf[i]->ppos[j-1]){
	fprintf(stderr,"\t-> Potential big problem, ordering of positions indicates that sites has not been sorted ? \n");
	fprintf(stderr,"\t-> chromoname: \'%s\' pos[%lu]:%d vs pos[%lu-1]:%d\n",chooseChr,j,saf[i]->ppos[j],j,saf[i]->ppos[j-1]);
	fprintf(stderr,"\t-> Consider running ./realSFS check your.saf.idx for each pop\n");
      }
      //      fprintf(stderr,"fl:%p fl->keeps:%p\n",fl,fl!=NULL?fl->keeps:NULL);
      if(i==0&&fl!=NULL&&fl->keeps!=NULL){//we only check for -sites for the first saf.
	
	if( saf[i]->toKeep->d[j] && fl->keeps[saf[i]->ppos[j]])
	  hit->d[saf[i]->ppos[j]]++;
      }else if(saf[i]->toKeep->d[j])
	hit->d[saf[i]->ppos[j]]++;
    }
  }

  if(killbreak){
    for(int i=0;i<saf.size();i++)
      saf[i]->dontRead =1;
    keep_destroy(hit);
    return 0;
  }
#if 0
  //  keep_info(hit,stderr,0,saf.size());
  for(int i=0;1&i<hit->m;i++)
    if(hit->d[i]==saf.size())
      fprintf(stdout,"%d\n",i+1);
  exit(0);
#endif
  //hit now contains the genomic position (that is the index).

  //let us now modify the the persaf toKeep char vector
  int tsk[saf.size()];
  int hasdata =0;
  for(int i=0;i<saf.size();i++) {
    tsk[i] =0;
    for(int j=0;j<=saf[i]->toKeep->last;j++)
      if(hit->d[saf[i]->ppos[j]]!=saf.size())
	saf[i]->toKeep->d[j] =0;
      else
	tsk[i]++;
    fprintf(stderr,"\t-> Sites to keep[%s] from pop%d:\t%d\n",chooseChr,i,tsk[i]);
    hasdata +=tsk[i];
    if(i>0)
      assert(tsk[i]==tsk[i-1]);
#if 0
    keep_info(saf[i]->toKeep,stderr,0,1);
    //print out overlapping positions for all pops
    for(int j=0;j<saf[i]->toKeep->last;j++){
      if(hit->d[saf[i]->ppos[j]]==saf.size())
	fprintf(stdout,"saf%d\t%d\n",i,j);
    }
#endif
  }
  if(hasdata==0){
    fprintf(stderr,"\t-> There is no data for this chr/scaffold lets skip\n");
    chooseChr=NULL;
    goto aGotoHereIsTheEasiest;
  }
  keep_destroy(hit);
  return 1;
}


size_t parspace(std::vector<persaf *> &saf){
  size_t ndim = 1;
  for(int i=0;i<saf.size();i++){
    ndim *= saf[i]->nChr+1;
    fprintf(stderr,"\t-> dim(%s):%lu\n",saf[i]->fname,saf[i]->nChr+1);
  }
  fprintf(stderr,"\t-> Dimension of parameter space: %lu\n",ndim);
  return ndim;
}



void setGloc(std::vector<persaf *> &saf,size_t nSites){
#if 1
  posiG = new int*[saf.size()];
  for(int i=0;i<saf.size();i++)
    posiG[i] = new int[nSites];
#endif
}

void delGloc(std::vector<persaf *> &saf,size_t nSites){
  for(int i=0;i<saf.size();i++)
    delete [] posiG[i];
  delete [] posiG;
}


size_t helper(persaf * pp,char *chr){
  if(chr==NULL)
    return pp->nSites;
  myMap::iterator it=pp->mm.find(chr);
  if(it==pp->mm.end()){
    fprintf(stderr,"\t-> Problem finding chromosome: %s\n",chr);
    exit(0);
  }
  return it->second.nSites;
}

/*
  returns the maxnumber of sites across all samples
 */
size_t calc_nsites(std::vector<persaf *> &pp,args *ar){
  if(ar->start!=-1 &&ar->stop!=-1)
    return ar->stop-ar->start;
  size_t res = helper(pp[0],ar->chooseChr);
  for(int i=1;i<pp.size();i++)
    if(helper(pp[i],ar->chooseChr) > res)
      res = helper(pp[i],ar->chooseChr);
  return res;
}




//just approximate
template <typename T>
size_t fsizes(std::vector<persaf *> &pp, size_t nSites){
  size_t res = 0;
  for(int i=0;i<pp.size();i++){
    res += nSites*(pp[i]->nChr+1)*sizeof(T)+nSites*sizeof( T*);
  }
  return res;
}

template size_t fsizes<float>(std::vector<persaf *> &pp, size_t nSites);


/*
  over elaborate function to read a sfs. Assumption is that input file contains the expected values.
  output is plugged into ret
 */

void readSFS(const char*fname,size_t hint,double *ret){
  fprintf(stderr,"\t-> Reading: %s assuming counts (will normalize to probs internally)\n",fname);
  FILE *fp = NULL;
  if(((fp=fopen(fname,"r")))==NULL){
    fprintf(stderr,"problems opening file:%s\n",fname);
    exit(0);
  }
  char buf[fsize(fname)+1];
  if(fsize(fname)!=fread(buf,sizeof(char),fsize(fname),fp)){
    fprintf(stderr,"Problems reading file: %s\n will exit\n",fname);
    exit(0);
  }
  buf[fsize(fname)]='\0';
  std::vector<double> res;
  char *tok=NULL;
  tok = strtok(buf,"\t\n ");
  if(!tok){
    fprintf(stderr,"File:%s looks empty\n",fname);
    exit(0);
  }
  res.push_back(atof(tok));

  while((tok=strtok(NULL,"\t\n "))) {  
    //fprintf(stderr,"%s\n",tok);
    res.push_back(atof(tok));

  }
  //  fprintf(stderr,"size of prior=%lu\n",res.size());
  if(hint!=res.size()){
    fprintf(stderr,"\t-> Pxroblem with size of dimension of prior %lu vs %lu\n",hint,res.size());
    for(size_t i=0;0&&i<res.size();i++)
      fprintf(stderr,"%zu=%f\n",i,res[i]);
    exit(0);
  }
  for(size_t i=0;i<res.size();i++){
      ret[i] = res[i];
      //      fprintf(stderr,"i=%lu %f\n",i,ret[i]);
  }
  normalize(ret,(int)res.size());
  for(int i=0;0&&i<res.size();i++)
    ret[i] = log(ret[i]);
  fclose(fp);
}


