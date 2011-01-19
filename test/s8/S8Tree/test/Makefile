CCC      = g++

# Subsystems that have compilable libraries
SUBSYS   = 

# Special rules to clean each subsystem
CLEANS   = $(addprefix CLN_,$(SUBSYS))

# List of programs with main functions to be filtered out of objects
PROGS    = main main_version

LIB		 = 

# Flags used in compilation
CXXFLAGS = -Wall -fPIC -I./ -I../ -I${ROOTSYS}/include -I${BOOST_ROOT}/include
LDFLAGS  = -shared -W1 `root-config --libs` -L../../lib

# Get list of all heads, sources and objects. Each source (%.cpp) whould have
# an object file except programs listed in PROGS
HEADS    = $(filter-out %_LinkDef.h,$(wildcard ./interface/*.h))
DICS	 = $(wildcard ./interface/*_LinkDef.h)
SRCS     = $(wildcard ./*.cpp)

OBJS     = $(foreach obj,$(addprefix ./obj/,$(patsubst %.cpp,%.o,$(notdir $(SRCS)))),$(filter-out $(addprefix %,$(addsuffix .o,$(PROGS))),$(obj)))
CINTS    = $(foreach dic,$(addprefix ./dic/,$(patsubst %.h,%.cxx,$(notdir $(DICS)))),$(subst _LinkDef,Dict,$(dic)))
CINTOBJS = $(foreach obj,$(subst ./dic/,,$(patsubst %.cxx,%.o,$(CINTS))),$(addprefix ./obj/,$(obj)))
LIBS     = -L../../lib $(addprefix -l,$(SUBSYS)) -lS8Tree

# Rules to be always executed: empty ones
.PHONY: $(PROGS)

all: $(SUBSYS) $(PROGS) lib

help:
	@echo "make <rule>"
	@echo
	@echo "Rules"
	@echo "-----"
	@echo
	@echo

lib: $(LIB)

$(PROGS): $(CINTS) $(CINTOBJS) $(OBJS) 
	@echo "[+] Compiling programs ..."
	$(CCC) $(CXXFLAGS) `root-config --glibs` -o ./$@ ./$@.cpp $(LIBS) $(OBJS) $(CINTOBJS)
	@echo

# Generate Dictionaries
$(CINTS): $(HEADS)
	@echo "[+] Generating ROOT Dictionaries ..."
	rootcint $@ -c $(HEADS) $(subst Dict.cxx,_LinkDef.h,$(subst dic,./interface,$@))
	@echo

$(CINTOBJS): $(CINTS)
	@echo "[+] Compiling ROOT Dictionaries ..."
	$(CCC) $(CXXFLAGS) -c $(addprefix ./dic/,$(patsubst %.o,%.cxx,$(notdir $@))) -o $@
	@echo

# Object files depend on all sources and headers but only sources should be
# compiled
$(OBJS): $(SRCS) $(HEADS)
	@echo "[+] Compiling objects ..."
	$(CCC) $(CXXFLAGS) -c $(addprefix ./src/,$(patsubst %.o,%.cpp,$(notdir $@))) -o $@
	@echo

$(LIB): $(CINTOBJS) $(OBJS)
	@echo "[+] Creating shared libraries ..."
	$(CCC) $(LDFLAGS) -o $@ $(OBJS) $(CINTOBJS) 
	@echo

# Call top level Makefile to compile subsystem this code depends on
$(SUBSYS):
	$(MAKE) -C ../ $@

# This rule will clean libraries also code depend on. Run:
#     make cleanall
# from top folder for global clean of all systems
cleanall: clean $(CLEANS)

$(CLEANS):
	$(MAKE) -C ../$(subst CLN_,,$@) clean

clean:
	rm -f ./obj/*.o
	rm -f ./dic/*.{h,cxx}
	rm -f ./lib/*.so
	rm -f $(addprefix ./,$(PROGS))
