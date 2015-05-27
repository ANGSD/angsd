CC=gcc
CXX=g++

FLAGS=-O3

HTSDIR = ../../htslib
HTSLIB = $(HTSDIR)/libhts.a



all: supersim thetaStat realSFS msToGlf smartCount printIcounts contamination splitgl NGSadmix

kmath.o: kmath.c kmath.h
	$(CC) $(FLAGS) -c kmath.c 

printIcounts: printIcounts.c
	$(CC) $(FLAGS) printIcounts.c -o printIcounts -lz

NGSadmix: ngsadmix32.cpp
	$(CXX) $(FLAGS) ngsadmix32.cpp -o NGSadmix -lz -lpthread

supersim: supersim.cpp rbeta.cpp
	$(CXX) $(FLAGS) supersim.cpp -lz -o supersim

contamination: contamination.cpp kmath.o
	$(CXX) $(FLAGS) contamination.cpp -lz -o contamination kmath.o -lpthread
splitgl: splitgl.c
	$(CC) $(FLAGS) splitgl.c -lz -o splitgl


safreader: safreader.cpp safreader.h
	$(CXX) $(FLAGS) -c safreader.cpp -lz -lpthread -I $(HTSDIR)

safstat: safstat.cpp safstat.h realSFS.h
	$(CXX) $(FLAGS) -c safstat.cpp -lz -lpthread -I $(HTSDIR)

realSFS_args: realSFS_args.cpp realSFS_args.h
	$(CXX) $(FLAGS) -c realSFS_args.cpp -I $(HTSDIR)


realSFS: realSFS.cpp safreader keep.hpp safstat realSFS_args
	$(CXX) $(FLAGS) realSFS.cpp -o realSFS  safreader.o realSFS_args.o safstat.o -I $(HTSDIR) $(HTSLIB) -lz -lpthread
	cp realSFS emOptim2

smartCount:	smartCount.c  
	$(CC) $(FLAGS) smartCount.c -I$(HTSDIR)/htslib   -lz -o smartCount -std=gnu99 -lpthread $(HTSLIB)

thetaStat:	thetaStat.cpp stats.cpp 
	$(CXX) $(FLAGS) thetaStat.cpp  -I$(HTSDIR)  -lz -o thetaStat -lpthread $(HTSLIB)
msToGlf: msToGlf6.c
	$(CC) $(FLAGS) msToGlf6.c -O3 -o msToGlf -std=gnu99 -lz -lm

clean:
	rm  -f *.o supersim emOptim emOptim2 msToGlf thetaStat smartCount realSFS printIcounts contamination spitgl NGSadmix splitgl *~
