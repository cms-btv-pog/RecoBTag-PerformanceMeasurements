CCC      = g++

SUBSYS   =

LIB		 = libStat.so

# Flags used in compilation
ifeq ($(strip $(DEBUG)),)
	DEBUG = -O2
else
	DEBUG = -O0 -g
endif

CXXFLAGS = ${DEBUG} -pipe -Wall -fPIC -I../ -I${ROOTSYS}/include -I${BOOST_ROOT}/include 
LDFLAGS  = -shared -W1 `root-config --libs` -L${BOOST_ROOT}/lib -L../lib $(foreach subsys,$(SUBSYS),$(addprefix -l,$(subsys)))

# Get list of all heads, sources and objects. Each source (%.cpp) whould have
# an object file
HEADS    = $(wildcard ./interface/*.h)
SRCS     = $(wildcard ./src/*.cpp)

OBJS     = $(foreach obj,$(addprefix ./obj/,$(patsubst %.cpp,%.o,$(notdir $(SRCS)))),$(obj))
LIBS     = -L${BOOST_ROOT}/lib

# Rules to be always executed: empty ones
.PHONY: all

all: lib

help:
	@echo "make <rule>"
	@echo
	@echo "Rules"
	@echo "-----"
	@echo
	@echo "  obj        compile objects"
	@echo "  lib        compile shared library"
	@echo

lib: $(LIB)

obj: $(OBJS)

# Object files depend on all sources and headers but only sources should be
# compiled
$(OBJS): $(addprefix ./src/,$(patsubst %.o,%.cpp,$(notdir $@))) 
	@echo "[+] Compiling objects ..."
	$(CCC) $(CXXFLAGS) -c $(addprefix ./src/,$(patsubst %.o,%.cpp,$(notdir $@))) -o $@
	@echo

$(LIB): $(CINTOBJS) $(OBJS)
	@echo "[+] Creating shared libraries ..."
	$(CCC) $(LDFLAGS) -o $(addprefix ./lib/,$@) $(OBJS) $(CINTOBJS) 
	@echo

# This rule will clean libraries also code depend on. Run:
#     make cleanall
# from top folder for global clean of all systems
cleanall: clean

clean:
	rm -f ./obj/*.o
	rm -f $(addprefix ./lib/,$(LIB))
