
CXX = $(CHARMC)
CXXFLAGS = -g --std=c++11

INC      = -I ${PETSC_DIR}/include -I ${PETSC_DIR}/${PETSC_ARCH}/include
INC     += -I ${PETSC_DIR}/include/mpiuni
LDLIBS  += -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc
CLDLIBS += -module CommonLBs
TLDLIBS += -lboost_unit_test_framework

DEPS  = region.hpp solver.hpp utils.hpp interpolator.hpp
OBJS  = region.o utils.o solver.o interpolator.o
TOBJS = testsuite.o test_utils.o test_solver.o test_region.o test_interpolator.o
COBJS = swr.o

all: $(OBJS)

.EXPORT_ALL_VARIABLES:
test:
	$(MAKE) -C tests test

testsuite: $(TOBJS) $(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS) $(TLDLIBS)

$(OBJS): %.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(TOBJS): %.o: tests/%.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

test-swr: swr
	./charmrun ++local +p4 ./swr 3 2 10 5 +balancer GreedyCommLB

swr: $(COBJS) $(OBJS)
	$(CHARMC) -O3 -language charm++ -o swr $^ $(LDLIBS) $(CLDLIBS)

SWRCharm.decl.h: swr.ci
	$(CHARMC)  swr.ci

swr.o: swr.cpp SWRCharm.decl.h
	$(CHARMC) $(CXXFLAGS) $(INC) -c swr.cpp

tags:
	etags *.cpp tests/*.cpp

clean:
	rm -vf nul *.o *.a tests/*.o
	rm -rf tests/testsuite swr charmrun
	rm -rf *.decl.h *.def.h

# For emacs flymake mode
check-syntax:
	$(CXX) $(CXXFLAGS) -c -o nul -Wall $(CHK_SOURCES)


