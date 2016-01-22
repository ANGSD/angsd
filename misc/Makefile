CC=gcc
CXX=g++

FLAGS=-O3

ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
$(info $(HTS_LIBDIR))
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif


all: supersim thetaStat realSFS msToGlf smartCount printIcounts contamination splitgl NGSadmix contamination2

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

contamination2: contamination2.cpp kmath.o
	$(CXX) $(FLAGS) contamination2.cpp -lz -o contamination2 kmath.o -lpthread

splitgl: splitgl.c
	$(CC) $(FLAGS) splitgl.c -lz -o splitgl

ifdef HTSSRC
safreader: safreader.cpp safreader.h header.h
	$(CXX) $(FLAGS) -c safreader.cpp -I $(HTS_INCDIR)

fstreader: fstreader.cpp fstreader.h header.h
	$(CXX) $(FLAGS) -c fstreader.cpp -I $(HTS_INCDIR)

safstat: safstat.cpp safstat.h 
	$(CXX) $(FLAGS) -c safstat.cpp -I $(HTS_INCDIR)

realSFS_args: realSFS_args.cpp realSFS_args.h
	$(CXX) $(FLAGS) -c realSFS_args.cpp -I $(HTS_INCDIR)

safcat: safcat.cpp safcat.h
	$(CXX) $(FLAGS) -c safcat.cpp -I $(HTS_INCDIR)

header: header.cpp header.h
	$(CXX) $(FLAGS) -c header.cpp -I $(HTS_INCDIR)

smartCount:	smartCount.c  
	$(CC) $(FLAGS) smartCount.c -I$(HTS_INCDIR)/htslib   -lz -o smartCount -std=gnu99 -lpthread $(HTS_LIBDIR)

thetaStat:	thetaStat.cpp stats.cpp 
	$(CXX) $(FLAGS) thetaStat.cpp  -I$(HTS_INCDIR)  -lz -o thetaStat -lpthread $(HTS_LIBDIR)

msToGlf: msToGlf6.c
	$(CC) $(FLAGS) msToGlf6.c  -I$(HTS_INCDIR) -O3 -o msToGlf -std=gnu99 -lz -lm  $(HTS_LIBDIR)


realSFS: realSFS.cpp safreader keep.hpp safstat realSFS_args header fstreader safcat
	$(CXX) $(FLAGS) realSFS.cpp -o realSFS  safreader.o realSFS_args.o safstat.o fstreader.o header.o safcat.o ../prep_sites.o ../analysisFunction.o -I $(HTS_INCDIR) $(HTS_LIBDIR) -lz -lpthread
else
safreader: safreader.cpp safreader.h header.h
	$(CXX) $(FLAGS) -c safreader.cpp 

fstreader: fstreader.cpp fstreader.h header.h
	$(CXX) $(FLAGS) -c fstreader.cpp 

safstat: safstat.cpp safstat.h 
	$(CXX) $(FLAGS) -c safstat.cpp 

realSFS_args: realSFS_args.cpp realSFS_args.h
	$(CXX) $(FLAGS) -c realSFS_args.cpp 

safcat: safcat.cpp safcat.h
	$(CXX) $(FLAGS) -c safcat.cpp

header: header.cpp header.h
	$(CXX) $(FLAGS) -c header.cpp 

smartCount:	smartCount.c  
	$(CC) $(FLAGS) smartCount.c -lz -o smartCount -std=gnu99 -lpthread -lhts

thetaStat:	thetaStat.cpp stats.cpp 
	$(CXX) $(FLAGS) thetaStat.cpp  -lz -o thetaStat -lpthread -lhts

msToGlf: msToGlf6.c
	$(CC) $(FLAGS) msToGlf6.c -O3 -o msToGlf -std=gnu99 -lz -lm  -lhts


realSFS: realSFS.cpp safreader keep.hpp safstat realSFS_args header fstreader safcat
	$(CXX) $(FLAGS) realSFS.cpp -o realSFS  safreader.o realSFS_args.o safstat.o fstreader.o header.o safcat.o ../prep_sites.o ../analysisFunction.o -lhts -lz -lpthread

endif

clean:
	rm  -f *.o supersim emOptim emOptim2 msToGlf thetaStat smartCount realSFS printIcounts contamination spitgl NGSadmix splitgl *~
