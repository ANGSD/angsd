CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)


FLAGS=-O3

PRG=htshook angsd misc


# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
HTSLIB = $(HTSDIR)/libhts.a
BGZIP  = $(HTSDIR)/bgzip

all: $(PRG)

.PHONY: misc clean htshook

misc:
	make -C misc/

htshook: string_alloc.h
	cmp string_alloc.h $(HTSDIR)/cram/string_alloc.h || cp -f string_alloc.h $(HTSDIR)/cram/string_alloc.h 
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

angsd.static: $(OBJ)
	$(CXX) $(FLAGS)  -o angsd.static *.o -lz -lpthread --static $(HTSLIB)

clean:
	rm  -f *.o *.d angsd angsd.static *~
	make -C misc/ clean
#	make -C ../htslib clean
test:
	echo "Not implemented yet"
