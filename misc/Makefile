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


all: supersim thetaStat realSFS msToGlf smartCount printIcounts contamination splitgl NGSadmix contamination2 haploToPlink ngsPSMC ibs scounts msHOT2glf

kmath.o: kmath.c kmath.h
	$(CC) $(FLAGS) -c kmath.c 

printIcounts: printIcounts.c
	$(CC) $(FLAGS) printIcounts.c -o printIcounts -lz

NGSadmix: ngsadmix32.cpp
	$(CXX) $(FLAGS) ngsadmix32.cpp -o NGSadmix -lz -lpthread

ibs: ibs.cpp
	$(CXX) $(FLAGS) ibs.cpp -o ibs -lz

supersim: supersim.cpp rbeta.cpp
	$(CXX) $(FLAGS) supersim.cpp -lz -o supersim

contamination: contamination.cpp kmath.o
	$(CXX) $(FLAGS) contamination.cpp -lz -o contamination kmath.o -lpthread

contamination2: contamination2.cpp kmath.o
	$(CXX) $(FLAGS) contamination2.cpp -lz -o contamination2 kmath.o -lpthread

haploToPlink: haploToPlink.cpp
	$(CXX) $(FLAGS) haploToPlink.cpp -lz -o haploToPlink

splitgl: splitgl.c
	$(CC) $(FLAGS) splitgl.c -lz -o splitgl

scounts: scounts.cpp
	$(CXX) $(FLAGS) scounts.cpp -lz -o scounts



ifdef HTSSRC
safreader.o: safreader.cpp safreader.h header.h
	$(CXX) $(FLAGS) -c safreader.cpp -I $(HTS_INCDIR)

psmcreader.o: psmcreader.cpp psmcreader.h header.h
	$(CXX) $(FLAGS) -c psmcreader.cpp -I $(HTS_INCDIR)

fstreader.o: fstreader.cpp fstreader.h header.h
	$(CXX) $(FLAGS) -c fstreader.cpp -I $(HTS_INCDIR)

safstat.o: safstat.cpp safstat.h 
	$(CXX) $(FLAGS) -c safstat.cpp -I $(HTS_INCDIR)

realSFS_optim.o: realSFS_optim.cpp realSFS_optim.h  multisafreader.hpp Matrix.hpp realSFS_shared.o
	$(CXX) $(FLAGS) -c realSFS_optim.cpp -I $(HTS_INCDIR)

realSFS_dadi.o: realSFS_dadi.cpp realSFS_dadi.h  multisafreader.hpp Matrix.hpp
	$(CXX) $(FLAGS) -c realSFS_dadi.cpp -I $(HTS_INCDIR)

realSFS_shared.o: realSFS_shared.cpp realSFS_shared.h 
	$(CXX) $(FLAGS) -c realSFS_shared.cpp -I $(HTS_INCDIR)

realSFS_args.o: realSFS_args.cpp realSFS_args.h safreader.o prep_sites
	$(CXX) $(FLAGS) -c realSFS_args.cpp -I $(HTS_INCDIR)

msHOT2glf: msHOT2glf.c
	$(CC) $(FLAGS) msHOT2glf.c  -I$(HTS_INCDIR) -O3 -o msHOT2glf -std=gnu99 $(HTS_LIBDIR)  -lz -lm -lbz2 -llzma -lpthread -lcurl

safcat.o: safcat.cpp safcat.h
	$(CXX) $(FLAGS) -c safcat.cpp -I $(HTS_INCDIR)

header.o: header.cpp header.h
	$(CXX) $(FLAGS) -c header.cpp -I $(HTS_INCDIR)

smartCount:	smartCount.c  
	$(CC) $(FLAGS) smartCount.c -I$(HTS_INCDIR)/htslib -o smartCount -std=gnu99 $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread -lcurl

thetaStat:	thetaStat.cpp stats.cpp 
	$(CXX) $(FLAGS) thetaStat.cpp  -I$(HTS_INCDIR)  -o thetaStat -lpthread $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread -lcurl

msToGlf: msToGlf.c
	$(CC) $(FLAGS) msToGlf.c  -I$(HTS_INCDIR) -O3 -o msToGlf -std=gnu99 $(HTS_LIBDIR)  -lz -lm -lbz2 -llzma -lpthread -lcurl

fpsmc.o: fpsmc.cpp fpsmc.h psmcreader.h main_psmc.h
	$(CXX) $(FLAGS) -c fpsmc.cpp -I $(HTS_INCDIR)

hmm_psmc.o: hmm_psmc.cpp hmm_psmc.h compute.c
	$(CXX) $(FLAGS) -c hmm_psmc.cpp -I $(HTS_INCDIR)


main_psmc.o: main_psmc.cpp main_psmc.h fpsmc.h hmm_psmc.o
	$(CXX) $(FLAGS) -c main_psmc.cpp -I$(HTS_INCDIR)

ngsPSMC: ngsPSMC.cpp psmcreader.o header.o main_psmc.o fpsmc.o hmm_psmc.o bfgs
	$(CXX) $(FLAGS) ngsPSMC.cpp  -O3 -o ngsPSMC psmcreader.o header.o main_psmc.o hmm_psmc.o fpsmc.o ../bfgs.o  -I$(HTS_INCDIR) $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread -lcurl

prep_sites: 
	make -C ../ prep_sites.o HTSSRC=$(realpath $(HTSSRC))

bfgs:
	make -C ../ bfgs.o HTSSRC=$(realpath $(HTSSRC))

analysisFunction: 
	make -C ../ analysisFunction.o HTSSRC=$(realpath $(HTSSRC))

realSFS: realSFS.cpp safreader.o keep.hpp safstat.o realSFS_args.o header.o fstreader.o safcat.o multisafreader.hpp prep_sites analysisFunction realSFS_optim.o realSFS_shared.o realSFS_dadi.o
	$(CXX) $(FLAGS) realSFS.cpp -o realSFS  safreader.o realSFS_args.o safstat.o fstreader.o header.o safcat.o realSFS_optim.o realSFS_dadi.o realSFS_shared.o ../prep_sites.o ../analysisFunction.o -I $(HTS_INCDIR) $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread
else
safreader.o: safreader.cpp safreader.h header.h
	$(CXX) $(FLAGS) -c safreader.cpp 

fstreader.o: fstreader.cpp fstreader.h header.h
	$(CXX) $(FLAGS) -c fstreader.cpp 

safstat.o: safstat.cpp safstat.h 
	$(CXX) $(FLAGS) -c safstat.cpp 

psmcreader.o: psmcreader.cpp psmcreader.h header.h
	$(CXX) $(FLAGS) -c psmcreader.cpp 

realSFS_args.o: realSFS_args.cpp realSFS_args.h safreader.o prep_sites
	$(CXX) $(FLAGS) -c realSFS_args.cpp 
safcat.o: safcat.cpp safcat.h
	$(CXX) $(FLAGS) -c safcat.cpp

realSFS_optim.o: realSFS_optim.cpp realSFS_optim.h realSFS_args.o multisafreader.hpp Matrix.hpp realSFS_shared.o
	$(CXX) $(FLAGS) -c realSFS_optim.cpp 

realSFS_dadi.o: realSFS_dadi.cpp realSFS_dadi.h realSFS_args.o multisafreader.hpp Matrix.hpp
	$(CXX) $(FLAGS) -c realSFS_dadi.cpp 

realSFS_shared.o: realSFS_shared.cpp realSFS_shared.h 
	$(CXX) $(FLAGS) -c realSFS_shared.cpp 

header.o: header.cpp header.h
	$(CXX) $(FLAGS) -c header.cpp 

smartCount:	smartCount.c  
	$(CC) $(FLAGS) smartCount.c -lz -o smartCount -std=gnu99 -lpthread -lhts  

thetaStat:	thetaStat.cpp stats.cpp 
	$(CXX) $(FLAGS) thetaStat.cpp  -lz -o thetaStat -lpthread -lhts -lz

msToGlf: msToGlf.c
	$(CC) $(FLAGS) msToGlf.c -O3 -o msToGlf -std=gnu99 -lhts  -lz -lm -lbz2 -llzma -lpthread

msHOTGlf: msHOTglf.c
	$(CC) $(FLAGS) msHOT2glf.c -O3 -o msHOT2glf -std=gnu99 -lhts  -lz -lm -lbz2 -llzma -lpthread    

prep_sites: 
	make -C ../ prep_sites.o

bfgs:
	make -C ../ bfgs.o

analysisFunction: 
	make -C ../ analysisFunction.o 

multisafreader.hpp: Matrix.hpp realSFS_shared.o 


realSFS: realSFS.cpp safreader.o keep.hpp safstat.o realSFS_args.o header.o fstreader.o safcat.o multisafreader.hpp prep_sites analysisFunction realSFS_optim.o realSFS_shared.o realSFS_dadi.o
	$(CXX) $(FLAGS) realSFS.cpp -o realSFS  safreader.o realSFS_args.o safstat.o fstreader.o header.o safcat.o realSFS_optim.o realSFS_shared.o ../prep_sites.o ../analysisFunction.o realSFS_dadi.o -lhts -lz -lpthread

hmm_psmc.o: hmm_psmc.cpp hmm_psmc.h compute.c
	$(CXX) $(FLAGS) -c hmm_psmc.cpp

fpsmc.o: fpsmc.cpp fpsmc.h psmcreader.h main_psmc.h
	$(CXX) $(FLAGS) -c fpsmc.cpp

main_psmc.o: main_psmc.cpp main_psmc.h fpsmc.h fpsmc.o
	$(CXX) $(FLAGS) -c main_psmc.cpp

ngsPSMC: ngsPSMC.cpp psmcreader.o header.o main_psmc.o fpsmc.o  hmm_psmc.o bfgs
	$(CXX) $(FLAGS) ngsPSMC.cpp -o ngsPSMC psmcreader.o header.o fpsmc.o main_psmc.o hmm_psmc.o ../bfgs.o -lhts -lz

endif

clean:
	rm  -f *.o supersim emOptim emOptim2 msToGlf thetaStat smartCount realSFS printIcounts contamination spitgl NGSadmix splitgl scounts *~
