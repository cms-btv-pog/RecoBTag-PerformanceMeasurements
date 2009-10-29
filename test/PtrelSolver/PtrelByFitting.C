void PtrelByFitting()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by Counting
    PtrelByFitting solver("templates_InclusiveMu5Pt50_SSVL.root", Fit::histograms);

    // Choose which templates to use for fitting
    // solver.useThreeTemplates();

    // Uncomment when using tagged light templates
    // solver.setNoTaggedLightTemplate(true);

    // Measure the efficiencies
    solver.solve(
        "/uscms_data/d2/pratima/PM_results/ThreeBin/results_Pt50_SSVL.root",
        "fitting_InclusiveMu5Pt50_SSVL.root"
    );
}

