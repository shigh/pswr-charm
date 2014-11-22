
CXX = g++
CXXFLAGS  = -g --std=c++11

INC      = -I ${PETSC_DIR}/include -I ${PETSC_DIR}/${PETSC_ARCH}/include
LDLIBS  += -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc
TLDLIBS += -lboost_unit_test_framework

DEPS  = region.hpp solver.hpp utils.hpp
OBJS  = region.o utils.o solver.o
TOBJS = testsuite.o test_utils.o


all: $(OBJS)

test: testsuite
	./testsuite --log_level=test_suite

testsuite: $(TOBJS) $(OBJS)
	$(CXX)  -o $@ $^ $(LDLIBS) $(TLDLIBS)

$(OBJS): %.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(TOBJS): %.o: tests/%.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@


clean:
	rm -vf nul *.o *.a
	rm -rf testsuite

# For emacs flymake mode
check-syntax:
	$(CXX) $(CXXFLAGS) -c -o nul -Wall $(CHK_SOURCES)


