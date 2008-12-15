// PtrelSolver

#ifdef NOSCRAMV

#include "PtrelByCounting.h"
#include "PtrelBySystem4.h"
#else
#include "RecoBTag/PerformanceMeasurements/test/PtrelSolverByCounting/PtrelByCounting.h"
#include "RecoBTag/PerformanceMeasurements/test/PtrelSolverByCounting/PtrelBySystem4.h"
#endif


#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PtrelByCounting+;
#pragma link C++ class PtrelBySystem4+;

#endif
