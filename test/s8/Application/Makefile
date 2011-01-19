CCC      = g++

# Subsystems that have compilable libraries
SUBSYS   = Controller

# List of programs with main functions to be filtered out of objects
PROGS    = $(patsubst ./src/%.cpp,%,$(wildcard ./src/*.cpp))

# Flags used in compilation
ifeq ($(strip $(DEBUG)),)
	DEBUG = -O2
else
	DEBUG = -O0 -g
endif

CXXFLAGS = ${DEBUG} -pipe -Wall -I../ -I../FWCore -I${ROOTSYS}/include -I${BOOST_ROOT}/include

LIBS     = -L../lib -lFWCore -lController


# Rules to be always executed: empty ones
.PHONY: $(PROGS)

all: $(PROGS)

help:
	@echo "make <rule>"
	@echo
	@echo "Rules"
	@echo "-----"
	@echo
	@echo "  all        compile executables"
	@echo

$(PROGS): $(OBJS) 
	@echo "[+] Compiling programs ..."
	$(CCC) $(CXXFLAGS) `root-config --glibs` $(LIBS) ./src/$@.cpp -o ../bin/s8_$@
	@echo

clean:
ifeq ($(strip $(PROGS)),)
else
	rm -f $(addprefix ../bin/s8_,$(PROGS))
endif
