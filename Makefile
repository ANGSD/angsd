#modied from htslib makefile
FLAGS=-O3 -D__WITH_POOL__

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: htshook angsd misc


BAMDIR=""
BDIR=$(realpath $(BAMDIR))
# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info "HTSSRC defined")
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
$(info $(HTS_LIBDIR))
else
$(info "HTSSRC not defined, assuming systemwide installation -lhts")
endif

PACKAGE_VERSION  = 0.902

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

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

angsd: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o angsd *.o -lz -lpthread $(HTS_LIBDIR)
else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

angsd: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o angsd *.o -lz -lpthread -lhts
endif

testclean:
	rm -rf test/sfstest/output test/tajima/output test/*.log version.h test/temp.txt

clean:	testclean
	rm  -f *.o *.d angsd angsd.static version.h *~
	make -C misc/ clean

test:
	echo "Only subset of analyses is being tested"
	cd test;./testAll.sh ../angsd $(BDIR)
force:
