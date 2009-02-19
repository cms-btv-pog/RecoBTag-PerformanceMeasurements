void MakeTemplatesSummer2008()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Define the template function
    TF1 bTemplate("btemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);

    // Set the initial parameters
    bTemplate.SetParameters(1.7256, -4.19201, .0000207, .849781, 1251050);

    // Define the template function
    TF1 clTemplate("cltemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);

    // Set the initial parameters
    clTemplate.SetParameters(1.7256, -4.19201, .0000207, .849781, 1251050);

    // Create a PtrelTemplateMaker
    PtrelTemplateMaker maker;

    // Set rebinning options
    // maker.rebin(Dependency::pT, 3);
    // maker.rebin(Dependency::eta, 3);
    maker.rebin(Dependency::pT, 9);
    maker.rebin(Dependency::eta, 9);
    maker.rebin(Dependency::ptrel, 2);

    // Setting up the template function forms
    maker.function(Flavor::b, bTemplate);
    maker.function(Flavor::cl, clTemplate);

    // Make templates
    maker.make("/uscms/home/pratima/PerformanceMeasurements/CMSSW_2_2_2/src/RecoBTag/PerformanceMeasurements/test/results_InclusiveMu5Pt50.root", "templates_InclusiveMu5Pt50.root");
}
