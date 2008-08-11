
// PtrelSolver
// Author: Victor E. Bazterra, UIC (2008)

#include "PtrelSolverDependencies.h"

ClassImp(Dependency)

const char * Dependency::Label[] = { "p_{T}^{rel} [GeV/c]", "p_{T} [GeV/c]", "#eta" };
const char * Dependency::Name[] = { "ptrel", "pT", "eta" };
const long   Dependency::Dimension = 3;

ClassImp(Flavor)

const char * Flavor::Label[] = {"b", "Non-b"};
const char * Flavor::Name[] = { "b", "cl" };
const long   Flavor::Dimension = 2;
