CCC      = g++

# Subsystems that have compilable libraries
SUBSYS   = Redirector

# Special rules to clean each subsystem
CLEANS   = $(addprefix CLN_,$(SUBSYS))

# List of programs with main functions to be filtered out of objects
PROGS    = test_multi_buffer

LIB		 =

# Flags used in compilation
CXXFLAGS = -I../
LDFLAGS  = -shared -fPIC

# Get list of all heads, sources and objects. Each source (%.cc) whould have
# an object file except programs listed in PROGS
HEADS    =
SRCS     = $(wildcard ./src/*.cc)

OBJS     =
LIBS     = -L../lib $(addprefix -l,$(SUBSYS))

# Rules to be always executed: empty ones
.PHONY: $(PROGS)

all: $(SUBSYS) $(PROGS)

help:
	@echo "make <rule>"
	@echo
	@echo "Rules"
	@echo "-----"
	@echo
	@echo "  progs      compile shared library. (to be used by ROOT)"
	@echo

progs: $(PROGS)

$(PROGS): $(SUBSYS)
	@echo "[+] Compiling programs ..."
	$(CCC) $(CXXFLAGS) -o ./bin/$@ ./src/$@.cc $(LIBS) $(CINTOBJS)
	@echo

# Call top level Makefile to compile subsystem this code depends on
$(SUBSYS):
	$(MAKE) -C ../

# This rule will clean libraries also code depend on. Run:
#     make cleanall
# from top folder for global clean of all systems
cleanall: clean $(CLEANS)

$(CLEANS):
	$(MAKE) -C ../$(subst CLN_,,$@) clean

clean:
	rm -f ./obj/*.o
	rm -f $(addprefix ./bin/,$(PROGS))
