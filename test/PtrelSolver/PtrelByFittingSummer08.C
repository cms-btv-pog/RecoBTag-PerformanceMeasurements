void PtrelByFittingSummer08()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by Counting
    PtrelByFitting solver("templates_InclusiveMu5Pt50.root",Fit::histograms);
    
    // Choose which templates to use for fitting
    solver.setFitFlavor(Flavor::b);
    //solver.setFitFlavor(Flavor::cl);
    solver.setFitFlavor(Flavor::c);
    solver.setFitFlavor(Flavor::l);

    // Measure the efficiencies
    solver.solve(
        "results_InclusiveMu5Pt50.root",
        "fitting_InclusiveMu5Pt50.root"
    );
}

