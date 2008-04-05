#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepLRObsPlots.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepLRValPlots.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepSolutionFilter.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtBTagAnalysis.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepSolutionFilter.h"
//#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepLRObsPlots.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepLRValPlots.h"
#include "RecoBTag/PerformanceMeasurements/interface/TtObservables.h"
#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistency.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(TtDilepLRValPlots);
DEFINE_ANOTHER_FWK_MODULE(TtDilepLRObsPlots);
DEFINE_ANOTHER_FWK_MODULE(TtDilepSolutionFilter);
DEFINE_ANOTHER_FWK_MODULE(TtBTagAnalysis);
DEFINE_ANOTHER_FWK_MODULE(TtSemilepSolutionFilter);
//DEFINE_ANOTHER_FWK_MODULE(TtSemilepLRObsPlots);
DEFINE_ANOTHER_FWK_MODULE(TtSemilepLRValPlots);
DEFINE_ANOTHER_FWK_MODULE(TtObservables);
DEFINE_ANOTHER_FWK_MODULE(TtTagConsistency);
