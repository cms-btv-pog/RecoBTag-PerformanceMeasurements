# make file to compile analysis macros

all:	
	@cd src && $(MAKE)
	@cd test/Plotters && $(MAKE)
	@cd test/PtrelSolver && $(MAKE)
	@cd test/S8Tools && $(MAKE)
	@cd test/S8Solver && $(MAKE)

clean:        
	@cd src && $(MAKE) clean
	@cd test/Plotters && $(MAKE) clean
	@cd test/PtrelSolver && $(MAKE)	clean
	@cd test/S8Tools && $(MAKE) clean
	@cd test/S8Solver && $(MAKE) clean


