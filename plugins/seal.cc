#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepLRObsPlots.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepLRValPlots.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepSolutionFilter.h"
#include "RecoBTag/PerformanceMeasurements/plugins/TtBTagAnalysis.h"
#include "RecoBTag/PerformanceMeasurements/interface/TtObservables.h"
#include "RecoBTag/PerformanceMeasurements/interface/TtSemiLeptonicTagCounting.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(TtDilepLRValPlots);
DEFINE_ANOTHER_FWK_MODULE(TtDilepLRObsPlots);
DEFINE_ANOTHER_FWK_MODULE(TtDilepSolutionFilter);
DEFINE_ANOTHER_FWK_MODULE(TtBTagAnalysis);
DEFINE_ANOTHER_FWK_MODULE(TtObservables);
DEFINE_ANOTHER_FWK_MODULE(TtSemiLeptonicTagCounting);
