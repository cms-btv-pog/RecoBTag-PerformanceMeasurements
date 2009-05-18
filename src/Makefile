# make file to compile analysis macros
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXext -lX11

CXX           = gcc

CXXFLAGS      = -fPIC -Wall -Wno-deprecated -m32 -DNOSCRAMV
LDFLAGS       = -fPIC -Wl -m32

LD_LIBRARY_PATH:= $(ROOTSYS)/lib/:$(LD_LIBRARY_PATH)

INCLUDEDIR    = -I$(ROOTSYS)/include -I../interface 

HEADERS       = BTagBaseTrackEvent.h BTagEvent.h BTagHistograms.h BTagLeptonEvent.h BTagTrackEvent.h CFractionFitter.h
OBJS          = BTagBaseTrackEvent.o BTagEvent.o BTagHistograms.o BTagLeptonEvent.o BTagTrackEvent.o BTagEventDict.o CFractionFitter.o

all:	$(OBJS)
	$(CXX) -shared -o libBTagEvent.so  $(LDFLAGS) $(OBJS) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR) 
%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDEDIR)

%Dict.cc: %LinkDef.h
	@echo "Generating dictionary $@"
	@rootcint -f $@ -c -DNOSCRAMV -I../interface $(HEADERS) $^

clean:        
	@rm -f $(OBJS) core *~ *.o *.so 

