# make file to compile analysis macros
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXext -lX11

CXX           = gcc

CXXFLAGS      = -fPIC -Wall -Wno-deprecated -m32
LDFLAGS       = -fPIC -Wl -m32

LD_LIBRARY_PATH:= $(ROOTSYS)/lib/:$(LD_LIBRARY_PATH)

INCLUDEDIR    = -I$(ROOTSYS)/include -I./ 

HEADERS       = PtrelByCounting.h PtrelByFitting.h PtrelSolver.h PtrelSolverDependencies.h PtrelTemplateMaker.h PtrelUtils.h 
OBJS          = PtrelByCounting.o PtrelByFitting.o PtrelSolver.o PtrelSolverDependencies.o PtrelTemplateMaker.o PtrelUtils.o PtrelSolverDict.o

all:	$(OBJS)
	$(CXX) -shared  -o libpTrelTools.so  $(LDFLAGS) $(OBJS) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR) 
%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR)

%Dict.cc: %LinkDef.h
	@echo "Generating dictionary $@"
	@rootcint -f $@ -c -DNOSCRAMV $(HEADERS) $^

clean:        
	@rm -f $(OBJS) core *~ *.o *.so PtrelSolverDict.*

