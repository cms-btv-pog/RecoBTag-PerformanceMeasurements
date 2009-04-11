void PtrelByFittingSummer08()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by Counting
    PtrelByFitting solver("rebin_templates_InclusiveMu5Pt50.root", Fit::histograms);

    // Choose which templates to use for fitting
    solver.useThreeTemplates();

    // Uncomment when using tagged light templates
    // solver.setNoTaggedLightTemplate(true);

    // Measure the efficiencies
    solver.solve(
        "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Summer08/Jan_20_2009/InclusiveMu5Pt50_Summer08_IDEAL_V9_v1_GEN-SIM-RECO_2/Results/results_InclusiveMu5Pt50.root",
        "rebin_fitting_InclusiveMu5Pt50.root"
    );
}

