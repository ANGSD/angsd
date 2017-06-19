


template<typename T>
void readGL(persaf *fp,size_t nSites,size_t dim,Matrix<T> *ret,int *pp, int scale2norm){
  // ret->x=nSites;
  ret->y=dim;
  size_t i;
  for(i=ret->x;SIG_COND&&i<nSites;i++){
    if(i>0 &&(i% howOften)==0  )
      fprintf(stderr,"\r\t-> Has read %fmio sites now at: %lu      ",howOften/1e6,i);
    //
    int pos;
    size_t bytes_read= iter_read(fp,ret->mat[i],sizeof(T)*dim,&pos);//bgzf_read(fp,ret->mat[i],sizeof(T)*dim);
    if(pp!=NULL)//setpos
      pp[i] =pos;//
    if(bytes_read!=0 && bytes_read<sizeof(T)*dim){
      fprintf(stderr,"\t-> Problem reading chunk from file, please check nChr is correct, will exit \n");
      exit(0);
    }
    if(bytes_read==0)
      break;

    for(size_t j=0;scale2norm&&j<dim;j++)
      ret->mat[i][j] = exp(ret->mat[i][j]);
  }
  ret->x=i;
  if(SIG_COND==0)
    exit(0);
}




//returns the number of sites read
template<typename T>
size_t readGLS(std::vector<persaf *> &adolf,size_t nSites,std::vector< Matrix<T> *> &ret,int **posi,int scale2norm){
  size_t pre=ret[0]->x;
  for(int i=0;i<adolf.size();i++){
    readGL(adolf[i],nSites,adolf[i]->nChr+1,ret[i],posi!=NULL?posi[i]:NULL,scale2norm);
    //fprintf(stderr,"adolf:%d\t%lu posi:%d tak:%d\n",i,ret[i]->x,posi[i][0],tak);
  }
     
  return ret[0]->x-pre;
}

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
  for(int i=0;i<saf.size();i++) {
    tsk[i] =0;
    for(int j=0;j<=saf[i]->toKeep->last;j++)
      if(hit->d[saf[i]->ppos[j]]!=saf.size())
	saf[i]->toKeep->d[j] =0;
      else
	tsk[i]++;
    fprintf(stderr,"\t-> Sites to keep[%s] from pop%d:\t%d\n",chooseChr,i,tsk[i]);
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
  keep_destroy(hit);
  return 1;
}



template <typename T>
int readdata(std::vector<persaf *> &saf,std::vector<Matrix<T> *> &gls,size_t nSites,char *chooseChr,int start,int stop, int *pp,char **curChr,filt *fl,int scale2norm){
  static size_t lastread=0;
  extern int ** posiG;
  //  fprintf(stderr,"[%s] nSites:%d lastread:%d\n",__FUNCTION__,nSites,lastread);
  if(lastread==0 ){
    fprintf(stderr,"\t-> Done reading data from chromosome will prepare next chromosome\n");
    int ret = set_intersect_pos(saf,chooseChr,start,stop,curChr,fl); 
    //    fprintf(stderr,"[%s] ret:%d\n",__FUNCTION__,ret);
    //    exit(0);
    if(ret==-3)
      return -3;
  }

  lastread=readGLS(saf,nSites,gls,posiG,scale2norm);
  if(lastread>0&&saf.size()>1)
    fprintf(stderr,"\t-> [%s] lastread:%lu posi:%d\n",__FUNCTION__,lastread,posiG[0][0]);
#if 1 //<- below con be removed when we believe all is working
  if(saf.size()>1&&lastread!=0)
    for(int i=1;i<saf.size();i++){
      fprintf(stderr,"\t-> Comparing positions: %d with 0 has:%lu\n",i,gls[0]->x);
      if(memcmp(posiG[0],posiG[i],gls[0]->x*sizeof(int))!=0){
	fprintf(stderr,"SAF file is out of sync contact developer\n");
	for(int s=0;s<gls[0]->x;s++){
	  //	  fprintf(stderr,"s:%d\n",s);
	  if(posiG[0][s]!=posiG[i][s]){
	    fprintf(stderr,"Mismatch at s:%d which is i=0 vs i=%d with pos1:%d pos2:%d\n",s,i,posiG[0][s],posiG[i][s]);
	    exit(0);
	  }
	}
      }
    }
  //  fprintf(stderr,"Done checking\n");
#endif
  if(lastread==0)
    fprintf(stderr,"\t-> Only read nSites: %lu will therefore prepare next chromosome (or exit)\n",gls[0]->x);
  //fprintf(stderr,"readdata lastread:%d\n\n",lastread);
  // exit(0);
  if(pp!=NULL)
    for(int i=0;i<gls[0]->x;i++)
      pp[i] = posiG[0][i];
  if(chooseChr!=NULL&&lastread==0){
    //fprintf(stderr,"return -2\n");
    return -2;
  }
  else if(chooseChr==NULL &&lastread==0 )
    return -2;
  else
    return 1;
}
