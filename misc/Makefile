CC  ?= gcc
CXX ?= g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl

ifdef HTSSRC
$(info HTSSRC defined)
CPPFLAGS += -I$(realpath $(HTSSRC))
LIBS := $(realpath $(HTSSRC))/libhts.a $(LIBS)
else
$(info HTSSRC not defined, assuming systemwide installation)
LIBS += -lhts
endif

FLAGS = $(CPPFLAGS) -O3 $(LDFLAGS)

.PHONY: all clean install

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin

INSTALL = install
INSTALL_DIR = $(INSTALL) -dm0755
INSTALL_PROGRAM = $(INSTALL) -Dm0755

PROGRAMS = thetaStat realSFS

PROGRAMS_MISC = supersim msToGlf smartCount printIcounts contamination splitgl NGSadmix contamination2 haploToPlink ngsPSMC ibs scounts msHOT2glf

all: $(PROGRAMS) $(PROGRAMS_MISC)

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

safreader.o: safreader.cpp safreader.h header.h
	$(CXX) $(FLAGS) -c safreader.cpp 

fstreader.o: fstreader.cpp fstreader.h header.h
	$(CXX) $(FLAGS) -c fstreader.cpp 

safstat.o: safstat.cpp safstat.h 
	$(CXX) $(FLAGS) -c safstat.cpp 

psmcreader.o: psmcreader.cpp psmcreader.h header.h
	$(CXX) $(FLAGS) -c psmcreader.cpp 

realSFS_args.o: realSFS_args.cpp realSFS_args.h safreader.o
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
	$(CC) $(FLAGS) smartCount.c -lz -o smartCount -std=gnu99 -lpthread $(LIBS)

thetaStat:	thetaStat.cpp stats.cpp 
	$(CXX) $(FLAGS) thetaStat.cpp -o thetaStat -lpthread $(LIBS) -lz

msToGlf: msToGlf.c
	$(CC) $(FLAGS) msToGlf.c -O3 -o msToGlf -std=gnu99 $(LIBS) -lz -lm -lbz2 -llzma -lpthread

msHOT2glf: msHOT2glf.c
	$(CC) $(FLAGS) msHOT2glf.c -O3 -o msHOT2glf -std=gnu99 $(LIBS)  -lz -lm -lbz2 -llzma -lpthread

multisafreader.hpp: Matrix.hpp realSFS_shared.o 


realSFS: realSFS.cpp safreader.o keep.hpp safstat.o realSFS_args.o header.o fstreader.o safcat.o multisafreader.hpp realSFS_optim.o realSFS_shared.o realSFS_dadi.o
	$(CXX) $(FLAGS) realSFS.cpp -o realSFS  safreader.o realSFS_args.o safstat.o fstreader.o header.o safcat.o realSFS_optim.o realSFS_shared.o ../prep_sites.o ../analysisFunction.o ../chisquare.o realSFS_dadi.o $(LIBS) -lz -lpthread

hmm_psmc.o: hmm_psmc.cpp hmm_psmc.h compute.c
	$(CXX) $(FLAGS) -c hmm_psmc.cpp

fpsmc.o: fpsmc.cpp fpsmc.h psmcreader.h main_psmc.h
	$(CXX) $(FLAGS) -c fpsmc.cpp

main_psmc.o: main_psmc.cpp main_psmc.h fpsmc.h fpsmc.o
	$(CXX) $(FLAGS) -c main_psmc.cpp

ngsPSMC: ngsPSMC.cpp psmcreader.o header.o main_psmc.o fpsmc.o  hmm_psmc.o
	$(CXX) $(FLAGS) ngsPSMC.cpp -o ngsPSMC psmcreader.o header.o fpsmc.o main_psmc.o hmm_psmc.o ../bfgs.o $(LIBS) -lz

clean:
	rm -f *.o $(PROGRAMS) $(PROGRAMS_MISC) *~

install: $(PROGRAMS)
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(PROGRAMS) $(DESTDIR)$(bindir)

install-misc: $(PROGRAMS_MISC)
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(PROGRAMS_MISC) $(DESTDIR)$(bindir)
