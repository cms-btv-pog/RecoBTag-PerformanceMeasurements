void PtrelByCountingCSA07()
{
	// Load the PtrelSolver library
	gSystem->Load("libPtrelSolver.so");
	
	// Create a PtrelSolver by counting
	PtrelByCounting solver("templates_qcd.root");
	
	// Measure the efficiencies
	solver.solve(
	  "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Mar-19-2008/results_QCD.root",
	  "templates_qcd.root",	  
	  "counting_qcd.root"
	);
}

