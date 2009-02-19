
#include "PtrelSolverDependencies.h"

ClassImp(Fit)

const char * Fit::Label[] = { "Functions(TMinuit)", "Histograms(TFractionFitter)" };
const char * Fit::Name[] = { "functions", "histograms" };
const long   Fit::Dimension = 2;

ClassImp(Dependency)

const char * Dependency::Label[] = { "p_{T}^{rel} [GeV/c]", "p_{T} [GeV/c]", "#eta" };
const char * Dependency::Name[] = { "ptrel", "pT", "eta" };
const long   Dependency::Dimension = 3;

ClassImp(Flavor)

const char * Flavor::Label[] = {"Non-b", "b"};
const char * Flavor::Name[] = { "cl", "b"};
const long   Flavor::Dimension = 2;
