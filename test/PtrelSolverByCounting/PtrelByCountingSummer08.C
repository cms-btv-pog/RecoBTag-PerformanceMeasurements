void PtrelByCountingSummer08()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolverByCounting.so");

    // Create a PtrelSolver by Counting
    PtrelByCounting solver("templates_InclusiveMu5Pt50.root", Fit::functions);

    // Measure the efficiencies
    solver.solve(
        "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Summer08/Jan_20_2009/QCDpt30_Summer08_IDEAL_V9_v4_GEN-SIM-RECO/Results/results_QCDpt30.root", 
        "templates_QCDpt30.root",    
        "counting_QCDpt30.root"
    );
}

