/*
  code based on bam_md.c from SAMtools
  
 */

#include <cassert>
#include <cmath>
#include <ctype.h>
#include "bams.h"
#include "kprobaln.h"


static inline int bam_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f' || x == 'F') return 4;
	else return 0;
}


#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += bam_aux_type2size(type); \
	} while(0)



uint8_t *bam_aux_get(const aRead &b, const char tag[2])
{
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = getAuxStart(&b);
	while (s < b.vDat + b.block_size) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return s;
		__skip_tag(s);
	}
	return 0;
}
char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};


int bam_cap_mapQ(aRead &b,char *ref,int thres){
  //  fprintf(stderr,"thres=%d mapq=%d\n",thres,b.mapQ);
  assert(ref!=NULL);
  char *seq =(char *) getSeq(&b), *qual =(char *) getQuals(&b);
  uint32_t *cigar = getCig(&b);
  //	bam1_core_t *c = &b->core;
  int i, x, y, mm, q, len, clip_l, clip_q;
  double t;
  if (thres < 0) thres = 40; // set the default
  mm = q = len = clip_l = clip_q = 0;
  for (i = y = 0, x = b.pos; i < b.nCig; ++i) {
    int j, l = cigar[i]>>4, op = cigar[i]&0xf;
    if (op == BAM_CMATCH) {
      for (j = 0; j < l; ++j) {
	int z = y + j;
	int c1 = bam1_seqi(seq, z);
	//	fprintf(stderr,"%c ",ref[x+j]);
	int c2 = bam_nt16_table[(int)ref[x+j]];
	if (ref[x+j] == 0) break; // out of boundary
	if (c2 != 15 && c1 != 15 && qual[z] >= 13) { // not ambiguous
	  ++len;
	  if (c1 && c1 != c2 && qual[z] >= 13) { // mismatch
	    ++mm;
	    q += qual[z] > 33? 33 : qual[z];
	  }
	}
      }
      if (j < l) break;
      x += l; y += l; len += l;
    } else if (op == BAM_CDEL) {
      for (j = 0; j < l; ++j)
	if (ref[x+j] == 0) break;
      if (j < l) break;
      x += l;
    } else if (op == BAM_CSOFT_CLIP) {
      for (j = 0; j < l; ++j) clip_q += qual[y+j];
      clip_l += l;
      y += l;
    } else if (op == BAM_CHARD_CLIP) {
      clip_q += 13 * l;
      clip_l += l;
    } else if (op == BAM_CINS) y += l;
    else if (op == BAM_CREF_SKIP) x += l;
  }
  for (i = 0, t = 1; i < mm; ++i)
    t *= (double)len / (i+1);
  t = q - 4.343 * log(t) + clip_q / 5.;
  if (t > thres) return -1;
  if (t < 0) t = 0;
  t = sqrt((thres - t) / thres) * thres;
  //	fprint

  //  fprintf(stderr, "%s %lf %d\n",b.vDat, t, q);
  return (int)(t + .499);
}

/*
  So we have 3 options,
  1) use exisiting
  2) do simple baq
  3) do extended baq

  //we assume, that the exisiting bq/zq are the simple baq
 */
int bam_prob_realn_core(aRead &b, const char *ref, int type)
{
  //fprintf(stderr,"[%s] b.vDat=%s\n",__FUNCTION__,b.vDat);
  int k, i, bw, x, y, yb, ye, xb, xe;
  
  //  fprintf(stderr,"apply_baq=%d extend_baq=%d\n",apply_baq,extend_baq);
  uint32_t *cigar = getCig(&b);
  kpa_par_t conf = kpa_par_def;
  uint8_t *bq = 0, *zq = 0;
  uint8_t  *qual = getQuals(&b);
  
  if(b.flag&BAM_FUNMAP||b.l_seq==0) return -1;

  if ((bq = bam_aux_get(b, "BQ")) != 0) {
    ++bq;
    //fprintf(stderr,"not here\n");
  }
  if ((zq = bam_aux_get(b, "ZQ")) != 0 && *zq == 'Z') {
    ++zq;
    //fprintf(stderr,"not here\n");
  }
  if (bq && type==2){
      //  fprintf(stderr,"is not here\n");
      //  bam_aux_del(b, bq-1);
      bq = 0;
    }
  if (bq && zq) { // remove the ZQ tag
    // fprintf(stderr,"is not here\n");
    // bam_aux_del(b, zq-1);
    zq = 0;
  }
  if ((bq || zq)&&type==1) {
    assert(bq!=zq);
    if (bq) { // then convert BQ to ZQ
      for (i = 0; i < b.l_seq; ++i)
	qual[i] = qual[i] + 64 < bq[i]? 0 : qual[i] - ((int)bq[i] - 64);
    } else if (zq ) { // then convert ZQ to BQ
      for (i = 0; i < b.l_seq; ++i)
	qual[i] += (int)zq[i] - 64;
    }else{
      fprintf(stderr,"Somestring strange \n");
      exit(0);
    }
    return 0;
  }


	// find the start and end of the alignment	
	x = b.pos, y = 0, yb = ye = xb = xe = -1;
	for (k = 0; k < b.nCig; ++k) {
		int op, l;
		op = cigar[k]&0xf; l = cigar[k]>>4;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			if (yb < 0) yb = y;
			if (xb < 0) xb = x;
			ye = y + l; xe = x + l;
			x += l; y += l;
		} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
		else if (op == BAM_CDEL) x += l;
		else if (op == BAM_CREF_SKIP) return -1; // do nothing if there is a reference skip
	}
	// set bandwidth and the start and the end
	bw = 7;
	if (abs((xe - xb) - (ye - yb)) > bw)
		bw = abs((xe - xb) - (ye - yb)) + 3;
	conf.bw = bw;
	xb -= yb + bw/2; if (xb < 0) xb = 0;
	xe += b.l_seq - ye + bw/2;
	if (xe - xb - b.l_seq > bw)
		xb += (xe - xb - b.l_seq - bw) / 2, xe -= (xe - xb - b.l_seq - bw) / 2;
	{ // glocal
		uint8_t *s, *r, *q,  *bq;
		uint8_t *seq = getSeq(&b);
		int *state;
		bq =(uint8_t*) calloc(b.l_seq + 1, 1);
		memcpy(bq, qual, b.l_seq);
		s =(uint8_t*) calloc(b.l_seq, 1);
		for (i = 0; i < b.l_seq; ++i) 
		       s[i] = bam_nt16_nt4_table[bam1_seqi(seq, i)];
		     r =(uint8_t*) calloc(xe - xb, 1);
		for (i = xb; i < xe; ++i) {
			if (ref[i] == 0) { xe = i; break; }
			r[i-xb] = bam_nt16_nt4_table[bam_nt16_table[(int)ref[i]]];
		}
		       state =(int *) calloc(b.l_seq, sizeof(int));
		     q =(uint8_t*) calloc(b.l_seq, 1);
		     kpa_glocal(r, xe-xb, s, b.l_seq, qual, &conf, state, q);
		if (type==1) { // in this block, bq[] is capped by base quality qual[]
		  //fprintf(stderr,"not extended\n");
			for (k = 0, x = b.pos, y = 0; k < b.nCig; ++k) {
				int op = cigar[k]&0xf, l = cigar[k]>>4;
				if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
					for (i = y; i < y + l; ++i) {
						if ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y)) bq[i] = 0;
						else bq[i] = bq[i] < q[i]? bq[i] : q[i];
					}
					x += l; y += l;
				} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
				else if (op == BAM_CDEL) x += l;
			}
			for (i = 0; i < b.l_seq; ++i) bq[i] = qual[i] - bq[i] + 64; // finalize BQ
		} else { // in this block, bq[] is BAQ that can be larger than qual[] (different from the above!)
		  //fprintf(stderr,"extended\n");
			uint8_t *left, *rght;
			left =(uint8_t*) calloc(b.l_seq, 1); rght =(uint8_t*) calloc(b.l_seq, 1);
			for (k = 0, x = b.pos, y = 0; k < b.nCig; ++k) {
				int op = cigar[k]&0xf, l = cigar[k]>>4;
				if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
					for (i = y; i < y + l; ++i)
						bq[i] = ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y))? 0 : q[i];
					for (left[y] = bq[y], i = y + 1; i < y + l; ++i)
						left[i] = bq[i] > left[i-1]? bq[i] : left[i-1];
					for (rght[y+l-1] = bq[y+l-1], i = y + l - 2; i >= y; --i)
						rght[i] = bq[i] > rght[i+1]? bq[i] : rght[i+1];
					for (i = y; i < y + l; ++i)
						bq[i] = left[i] < rght[i]? left[i] : rght[i];
					x += l; y += l;
				} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
				else if (op == BAM_CDEL) x += l;
			}
			for (i = 0; i < b.l_seq; ++i) bq[i] = 64 + (qual[i] <= bq[i]? 0 : qual[i] - bq[i]); // finalize BQ
			free(left); free(rght);
		}
		for (i = 0; i < b.l_seq; ++i) 
		  qual[i] -= bq[i] - 64; // modify qual
		free(bq); free(s); free(r); free(q); free(state);
	}
	return 0;
}
