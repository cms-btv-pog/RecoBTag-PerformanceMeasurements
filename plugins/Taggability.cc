//
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      Taggability
//
/**\class PerformanceMeasurements/Taggability

 Description:

	 Author: Francisco Yumiceva, Fermilab
*/
//
// $Id: Taggability.cc,v 1.5 2010/03/31 23:31:15 jindal Exp $
//
//

#include "RecoBTag/PerformanceMeasurements/interface/Taggability.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace edm;
using namespace std;

Taggability::Taggability(const edm::ParameterSet &iConfig)
{

    JetCollection_ = iConfig.getParameter<edm::InputTag>("JetCollection");
    useJetCorr_    = iConfig.getParameter<bool>("ApplyJetCorrections");
    jetCorrLabel_  = iConfig.getParameter<std::string>("jetCorrectionsLabel");
    MinJetPt_      = iConfig.getParameter<double>("MinPt");
    MaxJetEta_     = iConfig.getParameter<double>("MaxEta");
    MinNtrksInJet_ = iConfig.getParameter<int>("MinNtracks");
    MinTrkPtInJet_ = iConfig.getParameter<double>("MinTrkPt");
    MinNjets_      = iConfig.getParameter<int>("MinNjets");
    PVCollection_  = iConfig.getParameter<edm::InputTag>("PrimaryVertexCollection");
    MinNPV_        = iConfig.getParameter<int>("MinNPrimaryVertices");
    bTagTrackEventIPTagInfos_ = iConfig.getParameter<std::string>("bTagTrackEventIPtagInfos");
    writeHistos_   = iConfig.getParameter<bool>("WriteHistograms");

    edm::Service<TFileService> fs;
    h2_in  = fs->make<TH2F>("h2_in" ,"Jets without filter",10, MinJetPt_ , 150, 5, 0, MaxJetEta_ );
    h2_out = fs->make<TH2F>("h2_out","Taggability applied to jets",10, MinJetPt_ , 150, 5, 0, MaxJetEta_ );

}

Taggability::~Taggability() {}


bool Taggability::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    bool pass = false;

    // primary vertices
    Handle< View<reco::Vertex> > PVColl;
    iEvent.getByLabel(PVCollection_, PVColl);

    const View< reco::Vertex > &thePV = *PVColl;

    if (PVColl.isValid() == false )
    {
        edm::LogWarning("Taggability")
        <<" Some products not available in the event: VertexCollection "
        << PVCollection_<<" "
        << PVColl.isValid();
        return pass;
    }

    if ( (int)thePV.size() >= MinNPV_ ) pass = true;
    else return false;


    // Calo Jets
    Handle< View<reco::CaloJet> > jetsColl;
    iEvent.getByLabel(JetCollection_, jetsColl);

    const View< reco::CaloJet > &theJets = *jetsColl;

    // Get the bTagTrackEventIPTagInfo collection
    Handle<std::vector<reco::TrackIPTagInfo> > bTagTrackEventIPTagInfos;
    iEvent.getByLabel(bTagTrackEventIPTagInfos_, bTagTrackEventIPTagInfos);

    // initialize jet corrector
    const JetCorrector *acorrector = 0;
    if (useJetCorr_ == true)
    {
        acorrector = JetCorrector::getJetCorrector(jetCorrLabel_,iSetup);
    }

    int jetIndex = 0;
    int Njets = 0;

    double jetpt = 0;
    double jeteta= 0;

    for (edm::View<reco::CaloJet>::const_iterator jet = theJets.begin(); jet!=theJets.end(); ++jet)
    {

        // get jet corrections
        double jetcorrection = 1.;
        if (useJetCorr_ == true)
        {
            jetcorrection =  acorrector->correction(*jet);
        }
        jetpt = (jet->pt() * jetcorrection );
        jeteta = std::abs( jet->eta() );

        // Jet quality cuts
        if ( jetpt <= MinJetPt_ || jeteta >= MaxJetEta_ )
        {
            jetIndex++;
            continue;
        }

        h2_in->Fill( jetpt, jeteta );

        // Get a vector of reference to the selected tracks in each jet
        reco::TrackRefVector tracks( (*bTagTrackEventIPTagInfos)[jetIndex].selectedTracks() );

        int NgoodTrks = 0;

        for ( size_t index=0; index < tracks.size(); index++ )
        {
            reco::TrackRef track = tracks[index];

            if ( track->pt() > MinTrkPtInJet_ ) NgoodTrks++;

        }

        if (NgoodTrks < MinNtrksInJet_ )
        {
            jetIndex++;
            continue;
        }

        jetIndex++;
        Njets++; // good jets

        h2_out->Fill( jetpt, jeteta );
    }

    if ( Njets >= MinNjets_ ) pass = true;
    else pass = false;

    return pass;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Taggability);
