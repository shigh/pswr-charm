
CXX = g++
CXXFLAGS  = -g
DEPS = region.hpp
OBJS = region.o

all: 

# pswr: objs

%.o: %.cpp %.hpp

clean:
	rm -vf nul *.o *.a

# For emacs flymake mode
check-syntax:
	$(CXX) $(CXXFLAGS) -c -o nul -Wall $(CHK_SOURCES)

