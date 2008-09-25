void PtrelBySystem4()
{
    // Load the PtrelSolver library
    gSystem->Load("libPtrelSolver.so");

    // Create a PtrelSolver by counting
    PtrelBySystem4 solver("templates_ttbar_21X.root");

    // Measure the efficiencies
    solver.solve("results.root","system4_ttbar_21X.root");
}

