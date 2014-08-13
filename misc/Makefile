CC=gcc
CXX=g++

FLAGS=-O3 -D_USE_KNETFILE


all: supersim thetaStat realSFS msToGlf smartCount printIcounts contamination splitgl

bfgs.o: bfgs.c bfgs.h
	$(CC) $(FLAGS) -c bfgs.c -Wno-vexing-parse

printIcounts: printIcounts.c
	$(CC) $(FLAGS) printIcounts.c -o printIcounts -lz

supersim: supersim.cpp rbeta.cpp
	$(CXX) $(FLAGS) supersim.cpp -lz -o supersim

contamination: contamination.cpp 
	$(CXX) $(FLAGS) contamination.cpp -lz -o contamination
splitgl: splitgl.c
	$(CC) $(FLAGS) splitgl.c -lz -o splitgl

realSFS: realSFS.cpp realSFS_saf2.cpp
	$(CXX) $(FLAGS) realSFS.cpp bfgs.cpp -o realSFS -lz -lpthread -Wno-write-strings
	cp realSFS emOptim2
smartCount:	smartCount.c ../bgzf.o ../knetfile.o
	$(CC) $(FLAGS) smartCount.c -I../  ../bgzf.o -lz -o smartCount -std=gnu99 ../knetfile.o

thetaStat:	thetaStat.cpp stats.cpp ../kstring.o ../bgzf.o ../knetfile.o
	$(CXX) $(FLAGS) thetaStat.cpp -I../  ../kstring.o ../bgzf.o ../knetfile.o -lz -o thetaStat
msToGlf: msToGlf6.c
	$(CC) $(FLAGS) msToGlf6.c -O3 -o msToGlf -std=gnu99 -lz -lm

clean:
	rm  -f *.o supersim emOptim emOptim2 msToGlf thetaStat smartCount realSFS printIcounts *~
