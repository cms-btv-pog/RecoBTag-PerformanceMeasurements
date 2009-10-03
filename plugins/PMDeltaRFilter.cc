//
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      PMDeltaRFilter
//
/**\class PerformanceMeasurements/PMDeltaRFilter

 Description:

	 Author: Francisco Yumiceva, Fermilab
*/
//
// $Id: PMDeltaRFilter.cc,v 1.1 2009/09/22 02:45:37 yumiceva Exp $
//
//

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Handle.h"
#include "RecoBTag/PerformanceMeasurements/interface/PMDeltaRFilter.h"

using namespace edm;
using namespace std;

PMDeltaRFilter::PMDeltaRFilter(const edm::ParameterSet &iConfig)
{
    jets_ = iConfig.getParameter<edm::InputTag>("Jets");
    muons_ = iConfig.getParameter<edm::InputTag>("Muons");
    MaxDeltaR_ = iConfig.getParameter<double>("MaxDeltaR");

}

bool PMDeltaRFilter::filter(edm::Event& iEvent , const edm::EventSetup & iSetup )
{

    bool answer = false;

    Handle< View<pat::Jet> > jetHandle;
    iEvent.getByLabel(jets_, jetHandle);

    Handle< View<pat::Muon> > muonHandle;
    iEvent.getByLabel(muons_, muonHandle);

    for (View< pat::Muon >::const_iterator muonIter = muonHandle->begin();
            muonIter != muonHandle->end(); ++muonIter)
    {

        if (answer) break;

        for (View< pat::Jet >::const_iterator jetIter = jetHandle->begin();
                jetIter != jetHandle->end(); ++jetIter)
        {

            if ( DeltaR<reco::Candidate>()( *muonIter, *jetIter ) < MaxDeltaR_ )
            {
                answer = true;
                break;
            }

        }

    }
    return answer;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PMDeltaRFilter);
