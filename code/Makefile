
CXX = g++
CXXFLAGS  = -g --std=c++11
DEPS = region.hpp solver.hpp utils.hpp
OBJS = region.o utils.o solver.o 
INC  = -I ${PETSC_DIR}/include -I ${PETSC_DIR}/${PETSC_ARCH}/include

all: $(OBJS)

#pswr: objs

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

clean:
	rm -vf nul *.o *.a

# For emacs flymake mode
check-syntax:
	$(CXX) $(CXXFLAGS) -c -o nul -Wall $(CHK_SOURCES)


