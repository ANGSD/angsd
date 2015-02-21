CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)

PLATFORM=$(shell uname )
FLAGS=-O3 -D_USE_KNETFILE

ifeq ($(PLATFORM),Darwin)
	PRG=htshook angsd misc
else
	PRG=htshook angsd angsd.static misc
endif

all: $(PRG)

.PHONY: misc clean htshook

misc:
	make -C misc/

htshook:
	make -C ../htslib

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
HTSLIB = $(HTSDIR)/libhts.a
BGZIP  = $(HTSDIR)/bgzip
include $(HTSDIR)/htslib.mk


%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTSDIR)/htslib $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTSDIR)/htslib $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTSDIR)/htslib $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTSDIR)/htslib $*.cpp >$*.d


angsd: $(OBJ)
	$(CXX) $(FLAGS)  -o angsd *.o -lz -lpthread $(HTSLIB)

angsd.static: $(OBJ)
	$(CXX) $(FLAGS)  -o angsd.static *.o -lz -lpthread --static $(HTSLIB)

clean:
	rm  -f *.o *.d angsd angsd.static *~
	make -C misc/ clean
	make -C ../htslib clean
test:
	echo "Not implemented yet"
