# QDP++ config file

#CONFIG=/u/home/wsun/software/package-9-10-18-cpu/parscalar/install/qdp++-double/bin/qdp++-config
CONFIG=/u/home/wsun/software/package-9-10-18-cpu/parscalar/install/chroma-double/bin/chroma-config
#CONFIG=/u/home/wsun/software/package-10-5-18-gpu/jit-llvm-nvptx/install/sm_37_omp/qdp++-double/bin/qdp++-config
#CONFIG=/u/home/wsun/software/package-10-5-18-gpu/jit-llvm-nvptx/install/sm_37_omp/chroma-double/bin/chroma-config
CXX=$(shell $(CONFIG) --cxx)
QDP_CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I. -DQUIET
QDP_LDFLAGS=$(shell $(CONFIG) --ldflags)
QDP_LIBS=$(shell $(CONFIG) --libs)

SOURCE=$(wildcard *.cc)
OBJS=$(patsubst %.cc,%.o,$(SOURCE))

all:main

main:$(OBJS)
	$(CXX) -o $@ $(QDP_CXXFLAGS) $^ $(QDP_LDFLAGS) $(QDP_LIBS)

%.o:%.cc
	$(CXX) -c -o $@ $(QDP_CXXFLAGS) $< $(QDP_LDFLAGS) $(QDP_LIBS)
clean:
	rm -rf main *.o
