
CXX = g++
CXXFLAGS  = -g
DEPS = region.hpp solver.hpp utils.hpp
OBJS = region.o utils.o solver.o 

all: 

#pswr: objs

%.o: %.cpp %.hpp

clean:
	rm -vf nul *.o *.a

# For emacs flymake mode
check-syntax:
	$(CXX) $(CXXFLAGS) -c -o nul -Wall $(CHK_SOURCES)


