#include <cassert>
#include "prep_sites.h"

#define BIN ".bin"
#define IDX ".idx"
/*
  Dumps and reads in the following index format
  chrId: int
  offset: int64
  len: int
  hasMajMin: 

 */



//return one+two
char *append(const char *one,const char*two){
  char *ret = new char[strlen(one)+strlen(two)+1];
  strcpy(ret,one);
  strcpy(ret+strlen(ret),two);
  return ret;
}

filt *filt_read(const char *fname){
  fprintf(stderr,"\t-> [%s] Reading binary representation:%s\n",__FILE__,fname);
  filt *ret = new filt;
  ret->bg =NULL;
  ret->fp =NULL;
  ret->keeps = ret->major = ret->minor=NULL;
  ret->hasMajMin =0;
  ret->curLen =0;
  char *bin_name=append(fname,BIN);
  char *idx_name=append(fname,IDX);

  if(isNewer(fname,bin_name)||isNewer(fname,idx_name)){ 
    fprintf(stderr,"\t-> Potential problem: File: \'%s\' looks newer than files: \'%s\',\'%s\'\n\t-> Please delete %s/%s files and rerun.\n",fname,bin_name,idx_name,bin_name,idx_name);
    exit(0);
  }

  ret->fp= fopen(idx_name,"r");
  ret->bg=bgzf_open(bin_name,"r");


  while(1){
    int chrId;
    if(0==fread(&chrId,sizeof(int),1,ret->fp))
      break;
    asdf_dats tmp;
    if(1!=fread(&tmp.offs,sizeof(int64_t),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
    if(1!=fread(&tmp.len,sizeof(int),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
    if(1!=fread(&ret->hasMajMin,sizeof(int),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
    ret->offs[chrId]=tmp;
  }
  
#if 0
  std::map<int,asdf_dats>::const_iterator it;
  fprintf(stderr,"hasMajoMin: %d\n",ret->hasMajMin);
  for(it=ret->offs.begin();0&&it!=ret->offs.end();++it)
    fprintf(stderr,"id:%d offs:%lld len:%d\n",it->first,it->second.offs,it->second.len);
#endif
  fprintf(stderr,"\t-> [%s] nChr: %zu loaded from binary filter file\n",__FILE__,ret->offs.size());
  if(ret->hasMajMin==1)
    fprintf(stderr,"\t-> Filterfile contains major/minor information\n");
  
  delete [] bin_name;
  delete [] idx_name;
  return ret;
}





void filt_gen(const char *fname,const aMap * revMap,const aHead *hd,int isBed){
  fprintf(stderr,"\t-> Filterfile: %s supplied will generate binary representations... \n",fname);
  if(isBed)
    fprintf(stderr,"\t-> Assuming that filter file is a bedfile\n");

  
  std::map<const char*,int,ltstr>::const_iterator it;
#if 0
  for(it= revMap->begin();it!=revMap->end();++it)
    fprintf(stderr,"%s->%d->%d\n",it->first,it->second,hd->l_ref[it->second]);
#endif

  gzFile gz = Z_NULL;
  gz = gzopen(fname,"r");
  if(gz==Z_NULL){
    fprintf(stderr,"Problem opening file:%s\n",fname);
    exit(0);
  }

  char* outnames_bin = append(fname,BIN);
  char* outnames_idx = append(fname,IDX);
    
  const char *delims = "\t \n";
  BGZF *cfpD = bgzf_open(outnames_bin,"w9");
  FILE *fp =NULL;
  fp=fopen(outnames_idx,"w");
  
  std::map <int,char> mm;//simple structure to check that input has been sorted by chr/contig
  char *ary = NULL;
  char *major = NULL;
  char *minor = NULL;
  int last=-1;
  int nCols = -1;
  int hasMajMin=0;
  char buf[LENS];

  extern int SIG_COND;
 
  while(SIG_COND && gzgets(gz,buf,LENS)){
    if(buf[0]=='#')
      continue;
    char chr[LENS] ;int id=-1;
    char maj='N';
    char min='N';
    
    int nRead=0;
    int posS=-1;int posE=-1;
    char fourthCol[LENS];
  
    if(isBed==0)
      nRead=sscanf(buf,"%s\t%d\t%c\t%c\n",chr,&posS,&maj,&min);
    else
      nRead=sscanf(buf,"%s\t%d\t%d\t%s\n",chr,&posS,&posE,fourthCol);
   
    assert(posS>=0);

    if(isBed==0)
      posS--;
  
    if(nRead!=2&&nRead!=4){
      fprintf(stderr,"\t-> Filterfile must have either 2 columns or 4 columns\n");
      SIG_COND=0;goto cleanup;
      exit(0);
    }
    if(posS<0||((isBed)>0&&(posE<posS))){
      fprintf(stderr,"\t-> Problem with entry in filterfile: col2:%d col3:%d\n",posS,posE);
      fprintf(stderr,"\t-> Offending line looks like:\'%s\'\n",buf);
      SIG_COND=0;goto cleanup;
      exit(0);
    }
    
    assert(nRead!=0);
    if(nCols==-1){
      nCols=nRead;
      if(nCols==4&&isBed==0)
	hasMajMin = 1;
      if(hasMajMin==1)
	fprintf(stderr,"\t-> Filterfile contains major/minor information\n");
  
    }

#if 0
    fprintf(stderr,"nRead=%d: %s %d %c %c\n",nRead,chr,posS,maj,min);
    fprintf(stderr,"nRead=%d: %s %d %d %f\n",nRead,chr,posS,posE,score);
    exit(0);
#endif    
    it = revMap->find(chr);
    if(it==revMap->end()){
      fprintf(stderr,"chr: %s from filterfile: %s doesn't exist in index, will exit()\n",chr,fname);
      exit(0);
    }else
      id=it->second;
    
    
    /*
      1. if we have observed a chromo change then dump data
      2. realloc
      3. plug in values
      (need to check that we have alloced a data <=> last!=-1)
    */
    if(last!=-1 &&last!=id){
      //      fprintf(stderr,"writing index: last=%d id=%d\n",last,id);
      assert(ary!=NULL);
      //write data and index stuff
      int64_t retVal =bgzf_tell(cfpD);
      fwrite(&last,1,sizeof(int),fp);
      fwrite(&retVal,1,sizeof(int64_t),fp);
      fwrite(&hd->l_ref[last],1,sizeof(int),fp);
      fwrite(&hasMajMin,1,sizeof(int),fp);
      bgzf_write(cfpD,ary,hd->l_ref[last]);//write len of chr
      if(hasMajMin){
	bgzf_write(cfpD,major,hd->l_ref[last]);//write len of chr
	bgzf_write(cfpD,minor,hd->l_ref[last]);//write len of chr
      }
    }
    if(last!=id){
      fprintf(stderr,"\t-> Parsing chr:\'%s\'\n",chr);
      std::map<int,char>::iterator it=mm.find(id);
      if(it!=mm.end()){
	fprintf(stderr,"\t-> Potential problem: Filter file: \'%s\', doesn't look sorted by chr\n",fname);
	fprintf(stderr,"\t-> Please sort: \'%s\' by column1, and remove the temporary files: \'%s/%s\' before rerunning\n",fname,outnames_bin,outnames_idx);
	exit(0);
      }else
	mm[id]=1;
      last=id;
      ary=new char[hd->l_ref[last]];
      memset(ary,0,hd->l_ref[last]);
      //      fprintf(stderr,"chr: %s id=%d len=%d\n",chr,last,hd->l_ref[last]);
      if(hasMajMin){
	major=new char[hd->l_ref[last]];
	memset(major,0,hd->l_ref[last]);
	minor=new char[hd->l_ref[last]];
	memset(minor,0,hd->l_ref[last]);
      }
      
    }
    if(posS > hd->l_ref[id]){
      fprintf(stderr,"Position in filter file:%s is after end of chromosome? Will exit\n",fname);
      exit(0);
    }
    if(isBed>0){
      // assert(posS<posE);
      if(posS>posE){
	fprintf(stderr,"Problem parsing bedfile, end position looks before start position: %d vs %d\n",posS,posE);
	exit(0);
      }
      for(int ii=posS;ii<posE;ii++)
	ary[ii]=1;
    }else{
      ary[posS] = 1;
      if(hasMajMin){
	major[posS] = refToInt[maj];
	minor[posS] = refToInt[min];
      }
    }
  }
  if(last!=-1){
    //fprintf(stderr,"writing index\n");
    assert(ary!=NULL);
    //write data and index stuff
    int64_t retVal =bgzf_tell(cfpD);
    fwrite(&last,1,sizeof(int),fp);
    fwrite(&retVal,1,sizeof(int64_t),fp);
    fwrite(&hd->l_ref[last],1,sizeof(int),fp);
    fwrite(&hasMajMin,1,sizeof(int),fp);
    bgzf_write(cfpD,ary,hd->l_ref[last]);//write len of chr
    if(hasMajMin){
      bgzf_write(cfpD,major,hd->l_ref[last]);//write len of chr
      bgzf_write(cfpD,minor,hd->l_ref[last]);//write len of chr
    }
  }

  fprintf(stderr,"\t-> Filtering complete: Observed: %zu different chromosomes from file:%s\n",mm.size(),fname);
  mm.clear();
 cleanup:
  if(gz) gzclose(gz);
  if(fp) fclose(fp);
  if(cfpD) bgzf_close(cfpD);
  
  if(SIG_COND==0){
    fprintf(stderr,"\n\t-> CTRL+C was detected, we will therefore assume that the build of the binary filterfiles are incomplete.\n\t-> Will therefore delete: \'%s\',\'%s\'\n",outnames_bin,outnames_idx);
    unlink(outnames_bin);
    unlink(outnames_idx);
  }

}




filt *filt_init(const char *fname,const aMap* revMap,const aHead *hd,int isBed){
  char *bin_name=append(fname,BIN);
  char *idx_name=append(fname,IDX);
  //First


  if(!aio::fexists(bin_name)||!aio::fexists(idx_name))
    filt_gen(fname,revMap,hd,isBed);
  delete [] bin_name;
  delete [] idx_name;
  extern int SIG_COND;
  if(SIG_COND==0)
    return NULL;
  return  filt_read(fname);
}

