all:libpyhuddingec.so

CXX ?= g++
CXX_ := $(CXX)

BOOSTROOT=$(CONDA_PREFIX)

include MODER/Makefile

CXX := $(CXX_)
CXXFLAGS+= -fPIC


ALL_PAIRS_HUDDINGE_OBJS=MODER/common.o  MODER/probabilities.o MODER/parameters.o MODER/matrix_tools.o MODER/my_assert.o MODER/combinatorics.o MODER/multinomial_helper.o MODER/bndm.o MODER/orientation.o MODER/data.o MODER/iupac.o MODER/suffix_array_wrapper.o MODER/kmer_tools.o MODER/huddinge.o


libpyhuddingec.so: pyhuddinge.cpp $(ALL_PAIRS_HUDDINGE_OBJS)
	$(CXX) $(CXXFLAGS)  -shared -o libpyhuddingec.so pyhuddinge.cpp $(ALL_PAIRS_HUDDINGE_OBJS)
