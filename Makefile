CC  ?= gcc
CXX ?= g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl

CRYPTO_TRY=$(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto 2>/dev/null -o /dev/null; echo $$?)
ifeq "$(CRYPTO_TRY)" "0"
$(info Crypto library is available to link; adding -lcrypto to LIBS)
LIBS += -lcrypto
else
$(info Crypto library is not available to link; will not use -lcrypto)
endif


#if htslib source is defined
ifdef HTSSRC

#if hts source is set to systemwide
ifeq ($(HTSSRC),systemwide)
$(info HTSSRC set to systemwide; assuming systemwide installation)
LIBS += -lhts

else

#if hts source path is given
# Adjust $(HTSSRC) to point to your top-level htslib directory
$(info HTSSRC defined: $(HTSSRC))
CPPFLAGS += -I"$(realpath $(HTSSRC))"
LIBHTS := $(HTSSRC)/libhts.a
LIBS := $(LIBHTS) $(LIBS)

endif

#if htssrc not defined
else

$(info HTSSRC not defined; using ./htslib (newest htslib develop, fetched by make))
$(info Use `make HTSSRC=/path/to/htslib` to build angsd using a local htslib installation)
$(info Use `make HTSSRC=systemwide` to build angsd using the systemwide htslib installation)


HTSSRC := $(CURDIR)/htslib
CPPFLAGS += -I$(HTSSRC)
LIBHTS := $(HTSSRC)/libhts.a
LIBS := $(LIBHTS) $(LIBS)
HTSLIB_DEP := .htslib

endif

HTSLIB_URL    = https://github.com/samtools/htslib.git
HTSLIB_BRANCH = develop

.PHONY: .htslib

#always build against the newest htslib develop commit (full clone, so htslib reports its own version correctly):
# - no ./htslib: clone it
# - ./htslib is a git checkout: move it to the newest develop (keep the old one if offline)
# - ./htslib without .git (release tarball): build it as shipped
.htslib:
	@if [ ! -e htslib/Makefile ]; then \
	  git clone --branch $(HTSLIB_BRANCH) --recurse-submodules $(HTSLIB_URL) htslib || exit 1; \
	elif [ -e htslib/.git ]; then \
	  if git -C htslib fetch -q --tags $(HTSLIB_URL) $(HTSLIB_BRANCH); then \
	    git -C htslib checkout -q FETCH_HEAD && git -C htslib submodule update -q --init --recursive || exit 1; \
	  else \
	    echo "Warning: could not fetch htslib $(HTSLIB_BRANCH); building with the existing ./htslib"; \
	  fi; \
	fi
	$(MAKE) -C $(HTSSRC)



#modied from htslib makefile
FLAGS = -O3 -D__STDC_FORMAT_MACROS
CPPFLAGS := $(filter-out -DNDEBUG,$(CPPFLAGS))
FLAGS2 = $(CPPFLAGS) $(FLAGS) $(LDFLAGS)

CFLAGS := $(FLAGS2) $(CFLAGS)
CXXFLAGS := $(FLAGS2) $(CXXFLAGS)

#for compiling with ZSTD which is used for .bgen file format
ifeq ($(WITH_ZSTD),1)
LIBS += -lzstd
CXXFLAGS += -D__ZSTD__
endif

CSRC = $(wildcard *.c)
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin

INSTALL = install
INSTALL_DIR = $(INSTALL) -dm0755
INSTALL_PROGRAM = $(INSTALL) -m0755

#$(info CFLAGS=$(CFLAGS))
#$(info CXXFLAGS=$(CXXFLAGS))

PROGRAMS = angsd


all: $(PROGRAMS) misc

BAMDIR=""
BDIR=$(realpath $(BAMDIR))

#version, first one found of:
# 1) git checkout: git describe, e.g. 0.942 on a tag, 0.942-3-gabc1234(-dirty) after it
# 2) source tarball: the VERSION file, which git archive fills in (see .gitattributes)
# 3) unknown
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --tags --match '[0-9]*' --always --dirty 2>/dev/null)
endif
ifeq "$(PACKAGE_VERSION)" ""
PACKAGE_VERSION := $(shell sed -n 's/^\([0-9][^ ]*\).*/\1/p' VERSION 2>/dev/null)
endif
ifeq "$(PACKAGE_VERSION)" ""
PACKAGE_VERSION := unknown
endif
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))

version.h:
	echo '#define ANGSD_VERSION "$(PACKAGE_VERSION)"' > $@

.PHONY: all clean install install-all install-misc misc test

misc: analysisFunction.o bfgs.o prep_sites.o aio.o chisquare.o 
	$(MAKE) -C misc HTSSRC="$(realpath $(HTSSRC))"

-include $(OBJ:.o=.d)

#objects include htslib headers, so ./htslib must be in place first (also with make -j)
$(OBJ): | $(HTSLIB_DEP)

%.o: %.c
	$(CC) -c  $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d


#relink when htslib was rebuilt (e.g. make moved ./htslib to a newer develop commit)
ifneq "$(LIBHTS)" ""
$(LIBHTS): | $(HTSLIB_DEP) ;
endif

angsd: version.h $(OBJ) $(LIBHTS)
	$(CXX) $(FLAGS) $(LDFLAGS) -o angsd *.o $(LIBS)


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
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(PROGRAMS) $(DESTDIR)$(bindir)
	$(MAKE) -C misc HTSSRC="$(realpath $(HTSSRC))" install

install-misc: misc
	$(MAKE) -C misc HTSSRC="$(realpath $(HTSSRC))" install-misc

install-all: install install-misc
