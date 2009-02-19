void PtrelByCountingSummer08()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolverByCounting.so");

    // Create a PtrelSolver by Counting
    PtrelByCounting solver("templates_InclusiveMu5Pt50.root", Fit::histograms);

    // Measure the efficiencies
    solver.solve(
        "/uscms/home/pratima/PerformanceMeasurements/CMSSW_2_2_2/src/RecoBTag/PerformanceMeasurements/test/results_InclusiveMu5Pt50.root", 
        "templates_InclusiveMu5Pt50.root",    
        "counting_Inclusive.root"
    );
}

