CC  ?= gcc
CXX ?= g++

#modied from htslib makefile
FLAGS = $(CPPFLAGS) -O3 -D__WITH_POOL__ ${LDFLAGS}

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin

MKDIR_P = mkdir -p
INSTALL = install -p
INSTALL_DIR     = $(MKDIR_P) -m 755
INSTALL_PROGRAM = $(INSTALL)


PROGRAMS = angsd


all: angsd misc

BAMDIR=""
BDIR=$(realpath $(BAMDIR))
# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif

PACKAGE_VERSION  = 0.921

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define ANGSD_VERSION "$(PACKAGE_VERSION)"' > $@

.PHONY: misc clean test

misc:
	make -C misc/ HTSSRC=$(realpath $(HTSSRC))

-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c -I$(HTS_INCDIR)
	$(CC) -MM $(CFLAGS)  $*.c -I$(HTS_INCDIR)  >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR) 
	$(CXX) -MM $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR)  >$*.d

angsd: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o angsd *.o $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread -lcurl
else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

angsd: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o angsd *.o -lz -lpthread -lhts  -lbz2 -llzma
endif

testclean:
	rm -rf test/sfstest/output test/tajima/output test/*.log version.h test/temp.txt

clean:	testclean
	rm  -f *.o *.d angsd version.h *~
	make -C misc/ clean

test:
	echo "Only subset of analyses is being tested"
	cd test;./testAll.sh ../angsd $(BDIR)
force:

install: all
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(misc_bindir)
	$(INSTALL_PROGRAM) $(PROGRAMS) $(DESTDIR)$(bindir)
