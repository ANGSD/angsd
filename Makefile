CC  ?= gcc
CXX ?= g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
CPPFLAGS += -I$(realpath $(HTSSRC))
LIBS := $(realpath $(HTSSRC))/libhts.a $(LIBS)
else
$(info HTSSRC not defined, assuming systemwide installation)
LIBS += -lhts
endif

#modied from htslib makefile
FLAGS = $(CPPFLAGS) -O3 -D__WITH_POOL__ $(LDFLAGS)

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


all: $(PROGRAMS) misc

BAMDIR=""
BDIR=$(realpath $(BAMDIR))

PACKAGE_VERSION  = 0.921

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define ANGSD_VERSION "$(PACKAGE_VERSION)"' > $@

.PHONY: misc clean test

misc: analysisFunction.o bfgs.o prep_sites.o
	$(MAKE) -C misc HTSSRC=$(realpath $(HTSSRC))

-include $(OBJ:.o=.d)

%.o: %.c
	$(CC) -c  $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d

angsd: version.h $(OBJ)
	$(CXX) $(FLAGS) -o angsd *.o $(LIBS)

testclean:
	rm -rf test/sfstest/output test/tajima/output test/*.log version.h test/temp.txt

clean: testclean
	rm  -f *.o *.d $(PROGRAMS) version.h *~
	$(MAKE) -C misc clean

test:
	echo "Only subset of analyses is being tested"
	cd test;./testAll.sh ../angsd $(BDIR)
force:

install: all
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(misc_bindir)
	$(INSTALL_PROGRAM) $(PROGRAMS) $(DESTDIR)$(bindir)
