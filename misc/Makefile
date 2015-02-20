CC=gcc
CXX=g++

FLAGS=-O3 -D_USE_KNETFILE


all: supersim thetaStat realSFS msToGlf smartCount printIcounts contamination splitgl NGSadmix

bfgs.o: bfgs.c bfgs.h
	$(CC) $(FLAGS) -c bfgs.c -Wno-vexing-parse

printIcounts: printIcounts.c
	$(CC) $(FLAGS) printIcounts.c -o printIcounts -lz

NGSadmix: ngsadmix32.cpp
	$(CXX) $(FLAGS) ngsadmix32.cpp -o NGSadmix -lz -lpthread

supersim: supersim.cpp rbeta.cpp
	$(CXX) $(FLAGS) supersim.cpp -lz -o supersim

contamination: contamination.cpp 
	$(CXX) $(FLAGS) contamination.cpp -lz -o contamination
splitgl: splitgl.c
	$(CC) $(FLAGS) splitgl.c -lz -o splitgl

realSFS: realSFS.cpp realSFS_saf2.cpp
	$(CXX) $(FLAGS) realSFS.cpp bfgs.cpp -o realSFS -lz -lpthread -Wno-write-strings
	cp realSFS emOptim2
smartCount:	smartCount.c ../bgzf.o ../knetfile.o ../hfile.o ../hfile_net.o
	$(CC) $(FLAGS) smartCount.c -I../ ../hfile_net.o  ../bgzf.o -lz -o smartCount -std=gnu99 ../knetfile.o ../hfile.o -lpthread

thetaStat:	thetaStat.cpp stats.cpp ../kstring.o ../bgzf.o ../knetfile.o ../hfile.o ../hfile_net.o
	$(CXX) $(FLAGS) thetaStat.cpp -I../ ../hfile_net.o ../kstring.o ../bgzf.o ../knetfile.o ../hfile.o -lz -o thetaStat -lpthread
msToGlf: msToGlf6.c
	$(CC) $(FLAGS) msToGlf6.c -O3 -o msToGlf -std=gnu99 -lz -lm

clean:
	rm  -f *.o supersim emOptim emOptim2 msToGlf thetaStat smartCount realSFS printIcounts contamination spitgl NGSadmix splitgl *~
