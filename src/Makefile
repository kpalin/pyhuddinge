all:libpyhuddingec.so

include MODER/Makefile

CXXFLAGS+= -fPIC

ALL_PAIRS_HUDDINGE_OBJS=MODER/common.o  MODER/probabilities.o MODER/parameters.o MODER/matrix_tools.o MODER/my_assert.o MODER/combinatorics.o MODER/multinomial_helper.o MODER/bndm.o MODER/orientation.o MODER/data.o MODER/iupac.o MODER/suffix_array_wrapper.o MODER/kmer_tools.o MODER/huddinge.o
libpyhuddingec.so: pyhuddinge.cpp $(ALL_PAIRS_HUDDINGE_OBJS)
	g++ $(CXXFLAGS) $(LDFLAGS) -shared -o libpyhuddingec.so pyhuddinge.cpp $(ALL_PAIRS_HUDDINGE_OBJS)
