CCC      = g++

LIB		 = libRedirector.so

# Flags used in compilation
ifeq ($(strip $(DEBUG)),)
	DEBUG = -O2
else
	DEBUG = -O0 -g
endif

CXXFLAGS = ${DEBUG} -pipe -fPIC -Wall -I./ 
LDFLAGS  = -shared -W1

# Get list of all heads, sources and objects. Each source (%.cc) whould have
# an object file except programs listed in PROGS
HEADS    = $(wildcard ./interface/*.h)
SRCS     = $(wildcard ./src/*.cc)

OBJS     = $(foreach obj,$(addprefix ./obj/,$(patsubst %.cc,%.o,$(notdir $(SRCS)))),$(obj))
LIBS     =

# Rules to be always executed: empty ones
.PHONY: all

all: $(OBJS) lib

help:
	@echo "make <rule>"
	@echo
	@echo "Rules"
	@echo "-----"
	@echo
	@echo "  all        compile everythins."
	@echo "  lib        compile shared library."
	@echo

lib: $(LIB)

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

clean:
	@echo "[+] Cleaning ..."
	rm -f ./obj/*.o
	rm -f $(addprefix ./lib/,$(LIB))
	@echo
