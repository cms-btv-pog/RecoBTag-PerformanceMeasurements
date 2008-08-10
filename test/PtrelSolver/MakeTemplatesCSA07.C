void MakeTemplatesCSA07()
{
	// Load the PtrelSolver library
	gSystem->Load("libPtrelSolver.so");
	
	// Define the template function
	TF1 bTemplate("btemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);
	
	// Set the initial parameters
	bTemplate.SetParameters(1.30774, -0.51646, 0.00475143, 2.1, 800);

	// Define the template function
	TF1 clTemplate("cltemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);
	
	// Set the initial parameters
	clTemplate.SetParameters(1.30774, -2.51646, 0.00475143, 1.1, 800);

	// Create a PtrelTemplateMaker
	PtrelTemplateMaker maker;
	
	// Setting up the template function forms
	maker.function(Flavor::b, bTemplate);
	maker.function(Flavor::cl, clTemplate);
	
	// Make templates
	maker.make("/uscmst1b_scratch/lpc1/lpcbtag/pratima/Mar-19-2008/results_QCD.root", "templates_qcd.root");
}
