##modied from htslib makefile
CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)


FLAGS=-O3

all: htshook angsd misc


# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
HTS = $(realpath $(HTSDIR))
HTSLIB = $(HTS)/libhts.a

PACKAGE_VERSION  = 0.700

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define ANGSD_VERSION "$(PACKAGE_VERSION)"' > $@

.PHONY: misc clean htshook test

misc:
	make -C misc/ HTSDIR=$(HTS)

htshook: 
	make -C $(HTS)

-include $(OBJ:.o=.d)


%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS) $*.cpp >$*.d

angsd: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o angsd *.o -lz -lpthread $(HTSLIB)


testclean:
	rm -rf test/sfstest/output test/tajima/output test/*.log version.h

clean:	testclean
	rm  -f *.o *.d angsd angsd.static *~
	make -C misc/ clean

test:
	echo "Only subset of analyses is being tested"
	cd test;./testAll.sh ../angsd
force: