GECODE_LIBFLAGS=-lgecodekernel -lgecodesearch -lgecodeint -lgecodeminimodel  -lgecodesupport  -lgecodegist #-lGecodeExplorer

RNAAlignment_LIBS=$(GECODE_LIBFLAGS)

bin_PROGRAMS = RNAalignment exampleVisualization

exampleVisualization_SOURCES = exampleVisualization.cpp  \
	WinDisplay.cpp 

RNAalignment_SOURCES = RNAalignment.cpp alignmentScorePropagator/Neq.cc\
        alignmentScorePropagator/Neq.hh

