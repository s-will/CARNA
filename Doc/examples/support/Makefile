GECODE_HOME=/home/mmann/Software/gecode/3.0.2
GECODE_SRC_HOME=/home/mmann/Software/src/gecode-3.0.2

CPPFLAGS=-I$(GECODE_SRC_HOME)

CXX = g++
CXXFLAGS = $(CPPFLAGS) -Wextra -Wall -pipe -ggdb -fvisibility=hidden -O3 -fno-strict-aliasing -ffast-math $(CXXUSR)

LIBS=-L$(GECODE_HOME)/lib -lgecodekernel -lgecodesearch -lgecodeint -lgecodeminimodel #-lGecodeExplorer


all: example.o options.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@  $<

