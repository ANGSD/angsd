#pragma once


typedef struct  tNode_t{
  int l;//4
  int l2;//4//length of "insert"
  int m2;//4//possible length of insert before realloc
  int m;//4
  int refPos;//4
  char *seq;//8
  char *qs;//8
  unsigned char *posi;//8
  unsigned char *isop;//8
  unsigned char *mapQ;//8
  tNode_t *insert;//8 an insertion 8bytes;
  int deletion;//4 //counter 
}tNode;

typedef struct{
  int refId;
  int regStart;
  int regStop;
  int nSites;
  int nSamples;
  tNode ***nd;//nd[site][ind]
  int *refPos;//length is nSites
}chunkyT;



typedef struct{
  int l;//number of nodes in nodes
  int m; //possible number of of nodes
  int first;//simply a value which is equal to nodes[0].refPos;
  int last;
  tNode **nds;//this length can maximum be the maxlenght of a read.NOTANYMORE

}nodePoolT;


typedef struct{
  nodePoolT **dn;
  int nFiles;
  int regStart;
  int regStop;
  int refId;
}fcb;//forcallback

void dalloc_node(tNode *n);
chunkyT *mergeAllNodes_new(nodePoolT **dn,int nFiles);
void cleanUpChunkyT(chunkyT *chk);
