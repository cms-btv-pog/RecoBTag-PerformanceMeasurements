
#ifdef NOSCRAMV
#include "S8Solver.h"
#include "S8AnalyticSolver.h"
#else
#include "RecoBTag/PerformanceMeasurements/test/S8Solver/S8Solver.h"
#include "RecoBTag/PerformanceMeasurements/test/S8Solver/S8AnalyticSolver.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class S8Solver+;
#pragma link C++ class S8AnalyticSolver+;

#endif
