void Ptrel2DPlotter()
{
    // Load the Ptrel2DPlotter library
    gSystem->Load("libPlotters.so");

    // Declare the plotter
    Ptrel2DPlotter plotter("Ptrel2DPlotter.root");

    // Add a input root file
    plotter.Add("results.root");

    // Loop over the root
    plotter.Loop();

    plotter.Write();
}

