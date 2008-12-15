void PtrelByCountingCSA07_10invpb()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by counting
    PtrelByCounting solver("templates_qcd_WT_10invpb.root");

    // Measure the efficiencies
    solver.solve(
        "/uscmst1b_scratch/lpc1/lpcbtag/pratima/Mar-19-2008/results_QCD_WT_10pb.root",
        "templates_qcd_WT_10invpb.root",
        "counting_qcd_WT_10invpb.root"
    );
}

