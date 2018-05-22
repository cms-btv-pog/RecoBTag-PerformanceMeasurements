// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
// new includes
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class recoFilter : public edm::stream::EDFilter<> {
   public:
      explicit recoFilter(const edm::ParameterSet&);
   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      edm::EDGetTokenT<std::vector<pat::Muon>> muonsTok;
      edm::EDGetTokenT<std::vector<pat::Jet>> jetsAK8Tok;
      edm::EDGetTokenT<std::vector<pat::Jet>> jetsAK4Tok;
      edm::EDGetTokenT<std::vector<pat::Jet>> jetsAK4PuppiTok;
};

recoFilter::recoFilter(const edm::ParameterSet& iConfig)
{
   muonsTok = consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
   jetsAK8Tok = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJetsAK8"));
   jetsAK4Tok = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));
   jetsAK4PuppiTok = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJetsPuppi"));
}

bool recoFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<std::vector<pat::Muon>> muons;
   iEvent.getByToken(muonsTok, muons);
   bool nMuons = false;
   for (auto i = muons->begin(); i != muons->end(); ++i) {
      if (i->pt()>=50. && std::abs(i->eta())<2.1) {
         nMuons = true;
         break;
      }
   }
   if (!nMuons) return false;

   edm::Handle<std::vector<pat::Jet>> jetsAK8;
   iEvent.getByToken(jetsAK8Tok, jetsAK8);
   bool nJetsAK8 = false;
   for (auto i = jetsAK8->begin(); i != jetsAK8->end(); ++i) {
      if (i->pt()>=250. && std::abs(i->eta())<2.4) {
         nJetsAK8 = true;
         break;
      }
   }
   if (!nJetsAK8) return false;

   edm::Handle<std::vector<pat::Jet>> jetsAK4;
   iEvent.getByToken(jetsAK4Tok, jetsAK4);
   int nJetsAK4 = 0;
   for (auto i = jetsAK4->begin(); i != jetsAK4->end(); ++i) {
      if (i->pt()>=30. && std::abs(i->eta())<2.4) {
         ++nJetsAK4;
      }
   }
   
   edm::Handle<std::vector<pat::Jet>> jetsAK4Puppi;
   iEvent.getByToken(jetsAK4PuppiTok, jetsAK4Puppi);
   int nJetsAK4Puppi = 0;
   for (auto i = jetsAK4Puppi->begin(); i != jetsAK4Puppi->end(); ++i) {
      if (i->pt()>=30. && std::abs(i->eta())<2.4) {
         ++nJetsAK4Puppi;
      }
   }

   if (!(nJetsAK4>=2||nJetsAK4Puppi>=2)) return false;

   return true;
}

DEFINE_FWK_MODULE(recoFilter);

