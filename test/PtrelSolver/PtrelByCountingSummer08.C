void PtrelByCountingSummer08()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by Counting
    PtrelByCounting solver("templates_InclusiveMu5Pt50.root", Fit::histograms);

    // Measure the efficiencies
    solver.solve(
        "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Summer08/Jan_20_2009/InclusiveMu5Pt50_Summer08_IDEAL_V9_v1_GEN-SIM-RECO/Results/results_InclusiveMu5Pt50.root",
        "templates_InclusiveMu5Pt50.root",
        "counting_InclusiveMu5Pt50.root"
    );
}

