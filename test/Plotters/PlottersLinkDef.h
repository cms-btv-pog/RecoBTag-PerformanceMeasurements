#ifdef NOSCRAMV
#include "BasePlotter.h"
#include "Ptrel2DPlotter.h"
#else
#include "RecoBTag/PerformanceMeasurements/test/Plotters/BasePlotter.h"
#include "RecoBTag/PerformanceMeasurements/test/Plotters/Ptrel2DPlotter.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class BasePlotter+;
#pragma link C++ class Ptrel2DPlotter+;

#endif
