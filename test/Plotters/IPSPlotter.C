void IPSPlotter()
{
    gSystem->Load("libPlotters.so");
    
    IPSPlotter plot;
    plot.Add("results.root");
    plot.Loop(); // you can call this function like Loop(50000) limiting the number of events
    plot.Write();
}

