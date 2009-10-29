void PtrelTemplateMaker()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Define the template functions
    TF1 bTemplate("btemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);
    TF1 clTemplate("cltemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);
    TF1 cTemplate("ctemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);
    TF1 lTemplate("ltemplate", "[4] * (pow(x, [0]) * exp([1] * pow(x, [3])) + [2])", 0., 3.5);

    // Set the initial parameters
    bTemplate.SetParameters(1.7256, -4.19201, .0000207, .849781, 1251050);
    clTemplate.SetParameters(1.7256, -4.19201, .0000207, .849781, 1251050);
    cTemplate.SetParameters(1.7256, -4.19201, .0000207, .849781, 1251050);
    lTemplate.SetParameters(1.7256, -4.19201, .0000207, .849781, 1251050);

    // Create a PtrelTemplateMaker
    PtrelTemplateMaker maker;

    // Uncomment for running three bins approach
    // maker.rebin(Dependency::pT, 3);
    // maker.rebin(Dependency::eta, 3);

    // Uncomment for running one bin approach
    // maker.rebin(Dependency::pT, 9);
    // maker.rebin(Dependency::eta, 9);

    // Rebinning a factor to ptrel histograms
    maker.rebin(Dependency::ptrel, 2);

    // Setting up the template function forms
    maker.function(Flavor::b, bTemplate);
    maker.function(Flavor::cl, clTemplate);
    maker.function(Flavor::c, cTemplate);
    maker.function(Flavor::l, lTemplate);

    // Make templates
    maker.make(
        "/uscms_data/d2/pratima/PM_results/ThreeBin/results_Pt50_SSVL.root",
        "templates_InclusiveMu5Pt50_SSVL.root"
    );
}
