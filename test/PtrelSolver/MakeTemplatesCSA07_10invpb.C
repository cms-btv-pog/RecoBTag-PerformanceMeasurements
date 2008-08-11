void MakeTemplatesCSA07_10invpb()
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
    maker.rebin(Dependency::pT, 3);
    maker.rebin(Dependency::eta, 3);
    maker.rebin(Dependency::ptrel, 2);

    // Setting up the template function forms
    maker.function(Flavor::b, bTemplate);
    maker.function(Flavor::cl, clTemplate);

    // Make templates
    maker.make("/uscmst1b_scratch/lpc1/lpcbtag/pratima/Mar-19-2008/results_QCD_WT_10pb.root", "templates_qcd_WT_10invpb.root");
}
