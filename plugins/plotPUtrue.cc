// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
// new includes
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH1D.h>

class plotPUTrue : public edm::stream::EDAnalyzer<> {
   public:
      explicit plotPUTrue(const edm::ParameterSet&);
   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
      TH1D * h_getTrueNumInteractions;
};

plotPUTrue::plotPUTrue(const edm::ParameterSet& iConfig)
{
   puToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
   edm::Service<TFileService> fs;
   h_getTrueNumInteractions = fs->make<TH1D>("h_getTrueNumInteractions" , ";getTrueNumInteractions();events / 1", 100, 0., 100.);
}

void plotPUTrue::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<std::vector <PileupSummaryInfo>> puInfo;
   iEvent.getByToken(puToken_, puInfo);
   for (auto i = puInfo->begin(); i != puInfo->end(); ++i) {
      if (i->getBunchCrossing()==0) {
         h_getTrueNumInteractions->Fill(i->getTrueNumInteractions());
         break;
      }
   }
}

DEFINE_FWK_MODULE(plotPUTrue);

