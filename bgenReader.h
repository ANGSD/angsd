#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include "argStruct.h"
#include "analysisFunction.h"
#include "shared.h"

//reads .bgen files with layout 2 (https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html)
//only biallelic, unphased, diploid SNPs are used, other variants are skipped
class bgenReader{
private:
	FILE *fp;
	const char *fname;
	const aMap *revMap;
	int nInd;
	unsigned nVariants;//number of variants according to the header
	unsigned sitesRead;
	int compressed;//0: none, 1: zlib, 2: zstd

	std::vector<unsigned char> raw;//probability data as stored in the file
	std::vector<unsigned char> pds;//uncompressed probability data
	std::string chr;
	std::vector<std::string> alleles;

	//the current site
	int refId;
	int pos;//0-indexed
	char major;
	char minor;
	double *probs;//3*nInd

	int changed;//the current site is the first site of the next chromosome
	int nSkipWarn;
	int nPloidyWarn;

	void readHeader();
	void readBytes(void *dst,size_t size,const char *what);
	std::string readString(int lenBytes,const char *what);
	int readVariant();
	int decodeProbs();
	void copySite(funkyPars *r,int s);
public:
	funkyPars *fetch(int chunkSize);
	bgenReader(char *fname,const aMap *revMap_a,int intName_a,int &nInd_a);
	~bgenReader();
};
