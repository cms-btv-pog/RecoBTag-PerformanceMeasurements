CCC      = g++

LIB		 = libS8Tree.so

# Flags used in compilation
ifeq ($(strip $(DEBUG)),)
	DEBUG = -O2
else
	DEBUG = -O0 -g
endif

CXXFLAGS = ${DEBUG} -pipe -Wall -fPIC -I./ -I${ROOTSYS}/include -I${BOOST_ROOT}/include
LDFLAGS  = -shared -W1 `root-config --libs`

# Get list of all heads, sources and objects. Each source (%.cc) whould have
# an object file except programs listed in PROGS
HEADS    = $(wildcard ./interface/*.h)
DICS	 = $(wildcard ./src/*LinkDef.h)
SRCS     = $(wildcard ./src/*.cc)

OBJS     = $(foreach obj,$(addprefix ./obj/,$(patsubst %.cc,%.o,$(notdir $(SRCS)))),$(obj))
CINTS    = $(foreach dic,$(addprefix ./dic/,$(patsubst %.h,%.cxx,$(notdir $(DICS)))),$(subst LinkDef,Dict,$(dic)))
CINTOBJS = $(foreach obj,$(subst ./dic/,,$(patsubst %.cxx,%.o,$(CINTS))),$(addprefix ./obj/,$(obj)))

# Rules to be always executed: empty ones
.PHONY: $(PROGS)

all: dict $(OBJS) $(PROGS) lib

help:
	@echo "make <rule>"
	@echo
	@echo "Rules"
	@echo "-----"
	@echo
	@echo "  lib        compile shared library. (to be used by ROOT)"
	@echo

lib: $(LIB)

dict: $(CINTS) $(CINTOBJS)

# Generate Dictionaries
$(CINTS): $(HEADS)
	@echo "[+] Generating ROOT Dictionaries ..."
	rootcint $@ -c -I./ $(filter-out %Tools.h,$(HEADS)) $(subst Dict.cxx,LinkDef.h,$(subst dic,./src,$@))
	@echo

$(CINTOBJS): $(CINTS)
	@echo "[+] Compiling ROOT Dictionaries ..."
	$(CCC) $(CXXFLAGS) -c $(addprefix ./dic/,$(patsubst %.o,%.cxx,$(notdir $@))) -o $@
	@echo

# Object files depend on all sources and headers but only sources should be
# compiled
$(OBJS): $(SRCS) $(HEADS)
	@echo "[+] Compiling objects ..."
	$(CCC) $(CXXFLAGS) -c $(addprefix ./src/,$(patsubst %.o,%.cc,$(notdir $@))) -o $@
	@echo

$(LIB): $(CINTOBJS) $(OBJS)
	@echo "[+] Creating shared libraries ..."
	$(CCC) $(LDFLAGS) -o $(addprefix ./lib/,$@) $(OBJS) $(CINTOBJS) 
	@echo

# This rule will clean libraries also code depend on. Run:
#     make cleanall
# from top folder for global clean of all systems
cleanall: clean $(CLEANS)

$(CLEANS):
	$(MAKE) -C ../$(subst CLN_,,$@) clean

clean:
	rm -f ./obj/*.o
	rm -f ./dic/*.{h,cxx}
	rm -f $(addprefix ./lib/,$(LIB))
