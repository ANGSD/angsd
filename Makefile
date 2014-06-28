
FLAGS=-O3 -D_USE_KNETFILE


CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

all: angsd misc

.PHONY: misc clean

misc:
	make -C misc/

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)

%.o: %.c
	$(CC) -c  $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d


angsd: $(OBJ)
	$(CXX) $(FLAGS)  -o angsd *.o -lz -lpthread
	$(CXX) $(FLAGS)  -o angsd.static *.o -lz -lpthread --static

clean:
	rm  -f *.o *.d angsd angsd.static *~
	make -C misc/ clean
