void PtrelByFittingSummer08()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by Counting
    PtrelByFitting solver("onebin_templates_InclusiveMu5Pt50.root", Fit::histograms);

    // Choose which templates to use for fitting
    solver.setFitFlavor(Flavor::b);
    //solver.setFitFlavor(Flavor::cl);
    solver.setFitFlavor(Flavor::c);
    solver.setFitFlavor(Flavor::l);

    // Uncomment when using tagged light templates
    // solver.setTaggedLightTemplate(true);

    // Measure the efficiencies
    solver.solve(
        "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Summer08/Jan_20_2009/InclusiveMu5Pt50_Summer08_IDEAL_V9_v1_GEN-SIM-RECO_2/Results/results_InclusiveMu5Pt50.root",
        "onebin_fitting_InclusiveMu5Pt50.root"
    );
}

