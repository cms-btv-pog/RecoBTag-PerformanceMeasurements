#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepSolutionFilter.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtDilepSolutionFilter::TtDilepSolutionFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");
  weight = iConfig.getParameter< double > ("weight");
  jetPtCut   = iConfig.getParameter<double> ("jetPtCut");
  jetEtCut   = iConfig.getParameter<double> ("jetEtCut");
  leptonPtCut   = iConfig.getParameter<double> ("leptonPtCut");
  leptonTriggerPtCut   = iConfig.getParameter<double> ("leptonTriggerPtCut");

  evtsols    = iConfig.getParameter<edm::InputTag> ("EvtSolution");

  allSolution = 0;
  goodSolution = 0;
  B=0;nonB=0;tau=0;exception=0;
}


TtDilepSolutionFilter::~TtDilepSolutionFilter()
{
}

void TtDilepSolutionFilter::endJob()
{
   cout << "TtDilepSolutionFilter: All events " << allSolution <<" , "<< goodSolution <<" emu - "
   		   << B << " "<<nonB<<endl;

   cout << "TtDilepSolutionFilter: Events per fb-1 " << allSolution*weight <<" , "<< goodSolution*weight <<" emu - "
   		   << B*weight << " "<<nonB*weight<<" ; tau "<<tau*weight
		   <<" , "<< exception<<endl;
}

bool
TtDilepSolutionFilter::filter(edm::Event& iEvent, edm::EventSetup const & iSetup)
{
   // get the event solution
  edm::Handle< std::vector<TtDilepEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);
  std::vector<TtDilepEvtSolution> sols;
  sols = *eSols;

  if (debug) cout << "TtDilepSolutionFilter: Found "<< sols.size()  << " dilepton soutions.\n";

  //cout << "mass :";

  if (sols.size()!= 2)  return false;
  TtDilepEvtSolution & sol = sols[0];
 allSolution++;
  if ((sol.getCalJetB().pt()<jetPtCut) || (sol.getCalJetBbar().pt()<jetPtCut))
      return false;

  if ((sol.getCalJetB().et()<jetEtCut) || (sol.getCalJetBbar().et()<jetEtCut))
      return false;

  if (debug) cout << sol.getWpDecay() << " "<<sol.getWmDecay()<<endl;

  if ((sol.getLeptNeg().p4().pt()<leptonTriggerPtCut) && (sol.getLeptPos().p4().pt()<leptonTriggerPtCut))
	return false;

  if ((sol.getLeptNeg().p4().pt()<leptonPtCut) || (sol.getLeptPos().p4().pt()<leptonPtCut))
	return false;


  // Check that we have an e-mu event

  if ( ((sol.getWpDecay()!="muon")&&(sol.getWmDecay()!="muon")) ||
	((sol.getWpDecay()!="electron")&& (sol.getWmDecay()!="electron")) )
	return false;

//   if ( (sol.getWpDecay()=="electron")&&( !checkElectron(sol.getElectronp()) ) )
//      return false;
// 
//   if ( (sol.getWmDecay()=="electron")&&( !checkElectron(sol.getElectronm()) ) )
//      return false;

  ++goodSolution;

    try {if ( (abs(sols[0].getGenEvent()->lepton()->pdgId())==15) || 
  (abs(sols[0].getGenEvent()->leptonBar()->pdgId())==15) ) ++tau;
    } catch (...){exception++;}


    if (sols[0].getJetB().getPartonFlavour()==5) ++B;
      else ++nonB;
    if (sols[0].getJetBbar().getPartonFlavour()==5) ++B;
      else ++nonB;

 if (debug) cout <<" ============================ End TtDilepSolutionFilter ============================" <<endl;
  return true;
}


bool TtDilepSolutionFilter::checkElectron (const TopElectron & electron) const
{
  if ( (electron.hadronicOverEm() > 0.1) ||
  	(electron.eSuperClusterOverP() > 3.) ||
	(electron.eSuperClusterOverP() < 0.8) ) return false;
  return true;

}

