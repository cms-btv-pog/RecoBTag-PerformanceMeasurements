void PtrelBySystem4Summer08()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by Counting
    PtrelBySystem4 solver("rebin_templates_InclusiveMu5Pt50.root", Fit::histograms);

    // Uncomment when using three templates for fitting
    // solver.useThreeTemplates();
    
    // Uncomment when using untagged light templates
    // solver.setTaggedLightTemplate(false);

    // Measure the efficiencies
    solver.solve(
        "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Summer08/Jan_20_2009/InclusiveMu5Pt50_Summer08_IDEAL_V9_v1_GEN-SIM-RECO_2/Results/results_InclusiveMu5Pt50.root",
        "rebin_system4_InclusiveMu5Pt50.root"
    );
}

