// PtrelSolver

#ifdef NOSCRAMV

#include "PtrelByCounting.h"
#include "PtrelBySystem4.h"
#include "PtrelSolverDependencies.h"
#else
#include "RecoBTag/PerformanceMeasurements/test/PtrelSolverByCounting/PtrelByCounting.h"
#include "RecoBTag/PerformanceMeasurements/test/PtrelSolverByCounting/PtrelBySystem4.h"
#include "RecoBTag/PerformanceMeasurements/test/PtrelSolver/PtrelSolverDependencies.h"
#endif


#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Fit+;
#pragma link C++ class Flavor+;
#pragma link C++ class Dependency+;
#pragma link C++ class PtrelByCounting+;
#pragma link C++ class PtrelBySystem4+;

#endif
