#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <cstdarg>
#include <zlib.h>
#include "analysisFunction.h"
#include "bgenReader.h"

#ifdef __ZSTD__
#include <zstd.h>
#endif

static void bgenError(const char *fname,const char *fmt,...) __attribute__((format(printf,2,3)));
static void bgenError(const char *fname,const char *fmt,...){
	fprintf(stderr,"\t-> Problem with bgen file \'%s\': ",fname);
	va_list ap;
	va_start(ap,fmt);
	vfprintf(stderr,fmt,ap);
	va_end(ap);
	fprintf(stderr,". Will exit\n");
	exit(1);
}

void bgenReader::readBytes(void *dst,size_t size,const char *what){
	if(size>0&&fread(dst,1,size,fp)!=size)
		bgenError(fname,"file ends while reading %s of variant %u",what,sitesRead+1);
}

//string prefixed by its length stored in lenBytes (2 or 4) bytes
std::string bgenReader::readString(int lenBytes,const char *what){
	uint32_t len=0;
	if(lenBytes==2){
		uint16_t l16;
		readBytes(&l16,2,what);
		len=l16;
	}else
		readBytes(&len,4,what);
	std::string s(len,'\0');
	if(len>0)
		readBytes(&s[0],len,what);
	return s;
}

void bgenReader::readHeader(){
	uint32_t offset,LH,M,N,flags;
	if(fread(&offset,4,1,fp)!=1||fread(&LH,4,1,fp)!=1||fread(&M,4,1,fp)!=1||fread(&N,4,1,fp)!=1)
		bgenError(fname,"file is too short to contain a header");
	char magic[4];
	if(fread(magic,1,4,fp)!=4)
		bgenError(fname,"file is too short to contain a header");
	if(memcmp(magic,"bgen",4)&&memcmp(magic,"\0\0\0\0",4))
		bgenError(fname,"not a bgen file (magic number is not \'bgen\')");
	if(LH<20||LH>offset)
		bgenError(fname,"invalid header length %u (offset is %u)",LH,offset);
	if(fseek(fp,LH-20,SEEK_CUR)!=0||fread(&flags,4,1,fp)!=1)//skip the free data area
		bgenError(fname,"file is too short to contain a header");

	compressed = flags & 3;
	int layout = (flags>>2) & 15;
	if(layout!=2)
		bgenError(fname,"layout %d is not supported, only layout 2 is",layout);
	if(compressed==3)
		bgenError(fname,"unknown compression type 3");
#ifndef __ZSTD__
	if(compressed==2)
		bgenError(fname,"it is zstd compressed, but angsd is compiled without zstd (compile with WITH_ZSTD=1)");
#endif
	if(N==0||N>INT32_MAX/3)
		bgenError(fname,"invalid number of samples %u",N);
	nInd = N;
	nVariants = M;

	//sample identifiers (flags>>31) are not used, the variant data starts after offset+4 bytes
	if(fseek(fp,(long)offset+4,SEEK_SET)!=0)
		bgenError(fname,"cannot seek to the first variant");
}

bgenReader::bgenReader(char *fname_a,const aMap *revMap_a,int intName_a,int &nInd_a){
	fprintf(stderr,"\t-> NB -bgen input is experimental and only lightly tested. Only biallelic, unphased, diploid SNPs are used\n");
	fprintf(stderr,"\t-> NB Please report anything that looks strange at https://github.com/ANGSD/angsd/issues\n");
	fname = fname_a;
	revMap = revMap_a;
	fp = fopen(fname,"rb");
	if(fp==NULL){
		fprintf(stderr,"\t-> Problem opening bgen file: \'%s\'. Will exit\n",fname);
		exit(1);
	}
	readHeader();
	nInd_a = nInd;//number of samples comes from the file
	sitesRead = 0;
	probs = new double[3*nInd];
	refId = -1;
	pos = -1;
	major = minor = 4;
	changed = 0;
	nSkipWarn = 0;
	nPloidyWarn = 0;
}

bgenReader::~bgenReader(){
	fclose(fp);
	delete [] probs;
}

//decode the genotype probabilities in pds into probs.
//returns 1 if they are used, 0 if the site should be skipped
int bgenReader::decodeProbs(){
	size_t D = pds.size();
	if(D<10+(size_t)nInd)
		bgenError(fname,"genotype data of variant %u is too short",sitesRead);
	const unsigned char *p = &pds[0];
	uint32_t N2;memcpy(&N2,p,4);
	uint16_t K2;memcpy(&K2,p+4,2);
	if(N2!=(uint32_t)nInd)
		bgenError(fname,"variant %u has %u samples, the header says %d",sitesRead,N2,nInd);
	if(K2!=alleles.size())
		bgenError(fname,"variant %u has %u alleles in the genotype data and %zu in the variant data",sitesRead,K2,alleles.size());
	const unsigned char *pmis = p+8;//per sample: missing in the top bit, ploidy in the lower 6
	int phased = p[8+nInd];
	int B = p[9+nInd];
	if(B<1||B>32)
		bgenError(fname,"variant %u uses %d bits per probability, must be 1 to 32",sitesRead,B);
	if(phased){
		if(nSkipWarn++<10)
			fprintf(stderr,"\t-> Skipping phased site %s:%d (this message is printed at most 10 times)\n",chr.c_str(),pos+1);
		return 0;
	}

	//unphased biallelic: a sample with ploidy k has k stored values
	size_t nBits = 0;
	for(int i=0;i<nInd;i++)
		nBits += (size_t)(pmis[i]&63)*B;
	if(D<10+(size_t)nInd+(nBits+7)/8)
		bgenError(fname,"genotype data of variant %u is too short",sitesRead);

	const unsigned char *X = p+10+nInd;
	double scal = pow(2,B)-1;
	uint64_t mask = (B==32)?0xffffffffULL:((1ULL<<B)-1);
	size_t bit = 0;
	for(int i=0;i<nInd;i++){
		int ploidy = pmis[i]&63;
		int missing = pmis[i]>>7;
		uint32_t v[2]={0,0};
		for(int j=0;j<ploidy;j++){
			//read B bits starting at bit, they span at most 5 bytes
			uint64_t val=0;
			size_t byte=bit/8;
			size_t nBytes=(bit%8+B+7)/8;
			for(size_t b=0;b<nBytes;b++)
				val |= (uint64_t)X[byte+b]<<(8*b);
			if(j<2)
				v[j] = (uint32_t)((val>>(bit%8))&mask);
			bit += B;
		}
		double *pr = probs+3*i;
		if(missing||ploidy!=2){
			if(ploidy!=2&&!missing&&nPloidyWarn++<10)
				fprintf(stderr,"\t-> Sample %d has ploidy %d at %s:%d, it is treated as missing (this message is printed at most 10 times)\n",i,ploidy,chr.c_str(),pos+1);
			pr[0]=pr[1]=pr[2]=1.0/3.0;
			continue;
		}
		pr[0] = (1.0*v[0])/scal;
		pr[1] = (1.0*v[1])/scal;
		pr[2] = 1-(1.0*v[0])/scal-(1.0*v[1])/scal;
		if(pr[2]<0)//rounding
			pr[2] = 0;
		double sum = pr[0]+pr[1]+pr[2];
		pr[0] /= sum;
		pr[1] /= sum;
		pr[2] /= sum;
	}
	return 1;
}

//read the next variant. Returns 1 if it is a site that is used, 0 if it is skipped
int bgenReader::readVariant(){
	sitesRead++;
	readString(2,"the variant id");
	readString(2,"the rsid");
	chr = readString(2,"the chromosome");
	uint32_t vpos;
	uint16_t K;
	readBytes(&vpos,4,"the position");
	readBytes(&K,2,"the number of alleles");
	alleles.resize(K);
	for(int k=0;k<K;k++)
		alleles[k] = readString(4,"the alleles");

	uint32_t C,D;
	readBytes(&C,4,"the genotype data length");
	if(compressed){
		if(C<4)
			bgenError(fname,"invalid genotype data length %u for variant %u",C,sitesRead);
		readBytes(&D,4,"the genotype data length");
		C -= 4;
	}else
		D = C;
	raw.resize(C);
	readBytes(C?&raw[0]:NULL,C,"the genotype data");

	aMap::const_iterator it = revMap->find(chr.c_str());
	if(it==revMap->end())
		bgenError(fname,"chromosome \'%s\' of variant %u is not in the fai file",chr.c_str(),sitesRead);
	if(vpos<1||vpos>INT32_MAX)
		bgenError(fname,"invalid position %u for variant %u",vpos,sitesRead);
	refId = it->second;
	pos = vpos-1;

	if(K!=2||alleles[0].size()!=1||alleles[1].size()!=1){
		if(nSkipWarn++<10)
			fprintf(stderr,"\t-> Skipping site %s:%d, only biallelic SNPs are used (this message is printed at most 10 times)\n",chr.c_str(),pos+1);
		return 0;
	}
	major = refToChar[(unsigned char)alleles[0][0]];
	minor = refToChar[(unsigned char)alleles[1][0]];

	if(compressed==0)
		pds.swap(raw);
	else{
		pds.resize(D);
		if(compressed==1){
			uLongf destLen = D;
			int ret = uncompress(D?&pds[0]:NULL,&destLen,C?&raw[0]:NULL,C);
			if(ret!=Z_OK||destLen!=D)
				bgenError(fname,"cannot decompress the genotype data of variant %u (zlib error %d)",sitesRead,ret);
		}
#ifdef __ZSTD__
		else{
			size_t ret = ZSTD_decompress(D?&pds[0]:NULL,D,C?&raw[0]:NULL,C);
			if(ZSTD_isError(ret)||ret!=D)
				bgenError(fname,"cannot decompress the genotype data of variant %u (zstd)",sitesRead);
		}
#endif
	}
	return decodeProbs();
}

void bgenReader::copySite(funkyPars *r,int s){
	r->posi[s] = pos;
	r->major[s] = major;
	r->minor[s] = minor;
	r->keepSites[s] = nInd;
	memcpy(r->post[s],probs,3*nInd*sizeof(double));
}

funkyPars *bgenReader::fetch(int chunksize){
	funkyPars *r = funkyPars_init();
	r->posi = new int[chunksize];
	r->post = new double*[chunksize];
	r->major = new char[chunksize];
	r->minor = new char[chunksize];
	r->keepSites = new int[chunksize];
	r->nInd = nInd;
	for(int s=0;s<chunksize;s++)
		r->post[s] = new double[nInd*3];

	int nSites = 0;
	if(changed){//first site of this chromosome was read by the previous call
		r->refId = refId;
		copySite(r,nSites++);
		changed = 0;
	}
	while(nSites<chunksize&&sitesRead<nVariants){
		if(!readVariant())
			continue;
		if(nSites>0&&refId!=r->refId){
			changed = 1;
			break;
		}
		r->refId = refId;
		copySite(r,nSites++);
	}

	for(int s=nSites;s<chunksize;s++)
		delete [] r->post[s];
	r->numSites = nSites;
	if(nSites==0){
		funkyPars_destroy(r);
		return NULL;
	}
	return r;
}
