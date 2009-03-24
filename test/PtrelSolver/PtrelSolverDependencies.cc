
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

const char * Flavor::Label[] = { "b", "c", "l", "cl" };
const char * Flavor::Name[] = { "b", "c", "l", "cl" };
const long   Flavor::Dimension = 4;
