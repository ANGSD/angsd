
template<typename T>
void readGL(persaf *fp,size_t nSites,int dim,Matrix<T> *ret){

  // ret->x=nSites;
  ret->y=dim;
  size_t i;
  for(i=ret->x;SIG_COND&&i<nSites;i++){
    if(i>0 &&(i% howOften)==0  )
      fprintf(stderr,"\r\t-> Has read 5mio sites now at: %lu      ",i);
    //
    int pos;
    int bytes_read= iter_read(fp,ret->mat[i],sizeof(T)*dim,&pos);//bgzf_read(fp,ret->mat[i],sizeof(T)*dim);
    if(bytes_read!=0 && bytes_read<sizeof(T)*dim){
      fprintf(stderr,"Problem reading chunk from file, please check nChr is correct, will exit \n");
      exit(0);
    }
    if(bytes_read==0){
      //fprintf(stderr,"[%s] bytes_read==0 at i:%lu\n",__FUNCTION__,i);
      break;
    }
    for(size_t j=0;j<dim;j++)
      ret->mat[i][j] = exp(ret->mat[i][j]);
  }
  //fprintf(stderr,"[%s] i:%lu\n",__FUNCTION__,i);
  ret->x=i;
  if(SIG_COND==0)
    exit(0);
  //  matrix_print<T>(ret);
  //  exit(0);
  fprintf(stderr,"\r");
}



//returns the number of sites read
template<typename T>
int readGLS(std::vector<persaf *> &adolf,size_t nSites,std::vector< Matrix<T> *> &ret){
  int pre=ret[0]->x;
  for(int i=0;i<adolf.size();i++){
    readGL(adolf[i],nSites,adolf[i]->nChr+1,ret[i]);
    //    fprintf(stderr,"adolf:%d\t%lu\n",i,ret[i]->x);
  }
  
  return ret[0]->x-pre;
}

//unthreaded
//this will populate the keep vector by
// 1) set the chooseChr and populate toKeep
// 2) find over lap between different positions
// this is run once for each chromsome
int set_intersect_pos(std::vector<persaf *> &saf,char *chooseChr,int start,int stop){
  //  fprintf(stderr,"[%s] chooseChr:%s, start:%d stop:%d\n",__FUNCTION__,chooseChr,start,stop );

  if(saf.size()==1&&chooseChr==NULL){//use entire genome, then don't do any strange filtering
    //fprintf(stderr,"herer\n");
    return 0 ;
  }
  /*
    What happens here? very good question

   */
  static int firstTime =1;
  static myMap::iterator it_outer=saf[0]->mm.begin();
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
    fprintf(stderr,"\t-> Is in multi sfs, will now read data from chr:%s\n",chooseChr);
  }
  
  fprintf(stderr,"\t-> hello Im the master merge part of realSFS. and I'll now do a tripple bypass to find intersect \n");
  fprintf(stderr,"\t-> 1) Will set iter according to chooseChr and start and stop\n");
  assert(chooseChr!=NULL);

 //hit will contain the depth for the different populations
  keep<char> *hit =NULL;
  //  assert(saf.size()>1);

  if(saf.size()>1)
    hit =keep_alloc<char>();//
  
  
  //this loop will populate a 'hit' array containing the effective (differnt pops) depth
  //if we only have one population, then just return after iter_init
  for(int i=0;i<saf.size();i++){
    myMap::iterator it = iter_init(saf[i],chooseChr,start,stop);
    assert(it!=saf[i]->mm.end());  
    if(saf.size()==1)
      return 0;
    
    bgzf_seek(saf[i]->pos,it->second.pos,SEEK_SET);
    saf[i]->ppos = new int[it->second.nSites];
    bgzf_read(saf[i]->pos,saf[i]->ppos,it->second.nSites*sizeof(int));
    if(saf[i]->ppos[it->second.nSites-1] > hit->m)
      realloc(hit,saf[i]->ppos[it->second.nSites-1]+1);
    assert(hit->m>0);
    //    fprintf(stderr,"keep[%d].first:%lu last:%lu\n",i,saf[i]->toKeep->first,saf[i]->toKeep->last);
    for(int j=saf[i]->toKeep->first;j<saf[i]->toKeep->last;j++)
      if(saf[i]->toKeep->d[j])
	hit->d[saf[i]->ppos[j]]++;
  }
#if 0
  keep_info(hit,stderr,0,saf.size());
  for(int i=0;0&i<hit->m;i++)
    if(hit->d[i]==saf.size())
      fprintf(stdout,"%d\n",i);
  exit(0);
#endif
  //hit now contains the genomic position (that is the index).
  
  //let us now modify the the persaf toKeep char vector
  int tsk[saf.size()];
  for(int i=0;i<saf.size();i++){
    tsk[i] =0;
    for(int j=0;j<saf[i]->toKeep->last;j++)
      if(hit->d[saf[i]->ppos[j]]!=saf.size())
	saf[i]->toKeep->d[j] =0;
      else
	tsk[i]++;
    fprintf(stderr,"\t-> Sites to keep from pop%d:\t%d\n",i,tsk[i]);
    if(i>0)
      assert(tsk[i]==tsk[i-1]);
#if 0
    keep_info(saf[i]->toKeep,stderr,0,1);
    //print out overlapping posiitons for all pops
    
    for(int j=0;j<saf[i]->toKeep->last;j++){
      if(hit->d[saf[i]->ppos[j]]==saf.size())
	fprintf(stdout,"saf%d\t%d\n",i,j);
    }
#endif
  }
  //  exit(0);
  keep_destroy(hit);
}

template <typename T>
int readdata(std::vector<persaf *> &saf,std::vector<Matrix<T> *> &gls,int nSites,char *chooseChr,int start,int stop){

  static int lastread=0;
  //  fprintf(stderr,"[%s] nSites:%d lastread:%d\n",__FUNCTION__,nSites,lastread);
  if(lastread==0 ){
    //    fprintf(stderr,"\t-> Done reading data from chromosome will prepare next chromosome\n");
    int ret = set_intersect_pos(saf,chooseChr,start,stop); 
    //fprintf(stderr,"ret:%d\n",ret);
    if(ret==-3)
      return -3;
  }
  lastread=readGLS(saf,nSites,gls);
  if(lastread==0)
    fprintf(stderr,"\t-> Only read nSites: %lu will therefore prepare next chromosome (or exit)\n",gls[0]->x);
  //fprintf(stderr,"readdata lastread:%d\n\n",lastread);
  // exit(0);
  if(chooseChr!=NULL&&lastread==0){
    //fprintf(stderr,"return -2\n");
    return -2;
  }
  else if(chooseChr==NULL &&lastread==0 )
    return -2;
  else
    return 1;
}
