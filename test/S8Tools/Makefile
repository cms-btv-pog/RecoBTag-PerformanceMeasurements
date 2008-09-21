# make file to compile analysis macros
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

LIBS          = -L../../src -lBTagEvent $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXext -lX11

CXX           = gcc

CXXFLAGS      = -fPIC -Wall -Wno-deprecated -m32 -DNOSCRAMV
LDFLAGS       = -fPIC -Wl -m32

LD_LIBRARY_PATH:= $(ROOTSYS)/lib/:$(LD_LIBRARY_PATH)

INCLUDEDIR    = -I$(ROOTSYS)/include -I../../interface  -I./ 

HEADERS       = S8bPerformance.h  S8Plotter.h 
OBJS          = S8Plotter.o S8PlotterDict.o

all:	$(OBJS)
	$(CXX) -shared -o libS8Tools.so  $(LDFLAGS) $(OBJS) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR) 
%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR)

%Dict.cc: %LinkDef.h
	@echo "Generating dictionary $@"
	@rootcint -f $@ -c -I../../interface -DNOSCRAMV $(HEADERS) $^

clean:        
	@rm -f $(OBJS) core *~ *.o *.so PlottersDict.*

