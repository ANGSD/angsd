##modied from htslib makefile
CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)


FLAGS=-O3




# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
HTSLIB = $(HTSDIR)/libhts.a
BGZIP  = $(HTSDIR)/bgzip

PACKAGE_VERSION  = 0.700
NUMERIC_VERSION = $(PACKAGE_VERSION)

ifneq "$(wildcard .git)" ""
original_version := $(PACKAGE_VERSION)
PACKAGE_VERSION := $(shell git describe --always --dirty)
ifneq "$(subst ..,.,$(subst 0,,$(subst 1,,$(subst 2,,$(subst 3,,$(subst 4,,$(subst 5,,$(subst 6,,$(subst 7,,$(subst 8,,$(subst 9,,$(PACKAGE_VERSION))))))))))))" "."
empty :=
NUMERIC_VERSION := $(subst $(empty) ,.,$(wordlist 1,2,$(subst ., ,$(original_version))) 255)
endif

version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif



PRG=htshook angsd misc

all: $(PRG)



version.h:
	echo '#define ANGSD_VERSION "$(PACKAGE_VERSION)"' > $@

print-version:
	@echo $(PACKAGE_VERSION)

.PHONY: misc clean htshook test

misc:
	make -C misc/

htshook: 
	make -C $(HTSDIR)

-include $(OBJ:.o=.d)


%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTSDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTSDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTSDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTSDIR) $*.cpp >$*.d

angsd: $(OBJ)
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
