
CXX = $(CHARMC)
CXXFLAGS = -g

INC      = -I ${PETSC_DIR}/include -I ${PETSC_DIR}/${PETSC_ARCH}/include
INC     += -I ${PETSC_DIR}/include/mpiuni
LDLIBS  += -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc
CLDLIBS += -module CommonLBs
TLDLIBS += -lboost_unit_test_framework

DEPS   = region.hpp solver.hpp utils.hpp interpolator.hpp
OBJS   = region.o utils.o solver.o interpolator.o
TOBJS  = testsuite.o test_utils.o test_solver.o test_region.o test_interpolator.o
COBJS  = swr.o
CPOBJS = pswr.o

all: $(OBJS)

.EXPORT_ALL_VARIABLES:
test:
	$(MAKE) -C tests test

$(OBJS): %.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(TOBJS): %.o: tests/%.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

test-swr: swr
	./charmrun ++local +p4 ./swr 3 2 10 5 +balancer GreedyCommLB +LBDebug 3

test-pswr: pswr
	./charmrun ++local +p4 ./pswr 3 2 10 5 +balancer GreedyCommLB +LBDebug 3

swr: $(COBJS) $(OBJS)
	$(CHARMC) -O3 -language charm++ -o swr $^ $(LDLIBS) $(CLDLIBS) -tracemode projections

pswr: $(CPOBJS) $(OBJS)
	$(CHARMC) -O3 -language charm++ -o pswr $^ $(LDLIBS) $(CLDLIBS) -tracemode projections

SWRCharm.decl.h: swr.ci
	$(CHARMC)  swr.ci

swr.o: swr.cpp SWRCharm.decl.h
	$(CHARMC) $(CXXFLAGS) $(INC) -c swr.cpp

PSWRCharm.decl.h: pswr.ci
	$(CHARMC)  pswr.ci

pswr.o: pswr.cpp PSWRCharm.decl.h
	$(CHARMC) $(CXXFLAGS) $(INC) -c pswr.cpp

tags:
	etags *.cpp tests/*.cpp

clean:
	rm -vf nul *.o *.a tests/*.o
	rm -vf tests/testsuite swr pswr charmrun
	rm -vf *.decl.h *.def.h
	rm -vf *.log.gz *.projrc *.sts *.topo
	rm -vf ../taub/charmrun* ../taub/pswr_charm*

# For emacs flymake mode
check-syntax:
	$(CXX) $(CXXFLAGS) -c -o nul -Wall $(CHK_SOURCES)


