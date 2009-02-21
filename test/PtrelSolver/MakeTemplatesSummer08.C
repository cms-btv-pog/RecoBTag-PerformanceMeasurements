void MakeTemplatesSummer08()
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
    maker.rebin(Dependency::pT, 9);
    maker.rebin(Dependency::eta, 9);
    maker.rebin(Dependency::ptrel, 2);

    // Setting up the template function forms
    maker.function(Flavor::b, bTemplate);
    maker.function(Flavor::cl, clTemplate);

    // Make templates
    maker.make(
        "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Summer08/Jan_20_2009/QCDpt30_Summer08_IDEAL_V9_v4_GEN-SIM-RECO/Results/results_QCDpt30.root",
        "onebin_templates_QCDpt30.root"
    );
}
