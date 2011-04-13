# Subsystems
#
SUBSYS = Stat IO Utility Option Selector Analyzer Controller Application ROOT
SUBMOD = S8Tree Redirector FWCore

# Special rules to clean each subsystem
SS_CLEAN = $(addprefix CLN_,$(SUBSYS))
SM_CLEAN = $(addprefix CLN_,$(SUBMOD))

CLEANS = $(SS_CLEAN) $(SM_CLEAN)

# Rules to be always executed: empty ones
.PHONY: sys mod all $(SUBMOD) $(SUBSYS) cleanall $(CLEANS)

sys: $(SUBSYS)

mod: $(SUBMOD)

all: $(SUBMOD) $(SUBSYS)

help:
	@echo "make <rule>"
	@echo
	@echo "Rules"
	@echo "-----"
	@echo
	@echo "  all        Make all Subsystems: ${SUBSYS}."
	@echo "  clean      Soft clean (withouth submodules: ${SUBSYS}."
	@echo "  cleanall   Clean up all Subsystems including executables and libraries."
	@echo "  <subsys>   Make only one of the subsystems: ${SUBSYS}."
	@echo

$(SUBMOD):
	$(MAKE) -C $@ 
	ln -fs ../$@/lib/lib$@.so ./lib/

$(SUBSYS):
	$(MAKE) -C $@ 
	ln -fs ../$@/lib/lib$@.so ./lib/

clean: $(SS_CLEAN)

cleanall: $(CLEANS)

$(CLEANS):
	$(MAKE) -C $(subst CLN_,,$@) clean
	rm -f ./lib/lib$(subst CLN_,,$@).so
