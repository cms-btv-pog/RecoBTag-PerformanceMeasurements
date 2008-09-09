# make file to compile analysis macros
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

LIBS          = $(ROOTLIBS) -lMinuit2
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXext -lX11

CXX           = gcc

CXXFLAGS      = -fPIC -Wall -Wno-deprecated -m32
LDFLAGS       = -fPIC -Wl -m32

LD_LIBRARY_PATH:= $(ROOTSYS)/lib/:$(LD_LIBRARY_PATH)

INCLUDEDIR    = -I$(ROOTSYS)/include -I./ 

HEADERS       = S8AnalyticSolver.h S8fcn.h S8FitSolver.h S8NumericSolver.h S8Solver.h  
OBJS          = S8AnalyticSolver.o S8fcn.o S8FitSolver.o S8NumericSolver.o S8Solver.o S8SolverDict.o 

all:	$(OBJS)
	$(CXX) -shared  -o libS8Solver.so  $(LDFLAGS) $(OBJS) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR) 
%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR)

%Dict.cc: %LinkDef.h
	@echo "Generating dictionary $@"
	@rootcint -f $@ -c -DNOSCRAMV $(HEADERS) $^

clean:        
	@rm -f $(OBJS) core *~ *.o *.so PtrelSolverDict.*

