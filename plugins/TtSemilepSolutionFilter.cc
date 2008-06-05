#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepSolutionFilter.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtSemilepSolutionFilter::TtSemilepSolutionFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");
  weight = iConfig.getParameter< double > ("weight");
  jetPtCut   = iConfig.getParameter<double> ("jetPtCut");
  jetEtCut   = iConfig.getParameter<double> ("jetEtCut");
  leptonPtCut   = iConfig.getParameter<double> ("leptonPtCut");
  leptonTriggerPtCut   = iConfig.getParameter<double> ("leptonTriggerPtCut");
  leptonTrackIsoCut = iConfig.getParameter<double> ("leptonTrackIsoCut");
  leptonCaloIsoCut = iConfig.getParameter<double> ("leptonCaloIsoCut");


  evtsols    = iConfig.getParameter<edm::InputTag> ("EvtSolution");

  allSolution = 0;
  goodSolution = 0;
  B=0;nonB=0;tau=0;exception=0;
}


TtSemilepSolutionFilter::~TtSemilepSolutionFilter()
{
}

void TtSemilepSolutionFilter::endJob()
{
   cout << "TtSemilepSolutionFilter: All events " << allSolution <<" , "<< goodSolution <<" emu - "
   		   << B << " "<<nonB<<endl;

   cout << "TtSemilepSolutionFilter: Events per fb-1 " << allSolution*weight <<" , "<< goodSolution*weight <<" emu - "
   		   << B*weight << " "<<nonB*weight<<" ; tau "<<tau*weight
		   <<" , "<< exception<<endl;
}

bool
TtSemilepSolutionFilter::filter(edm::Event& iEvent, edm::EventSetup const & iSetup)
{
   // get the event solution
  edm::Handle< std::vector<TtSemiEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);
  std::vector<TtSemiEvtSolution> sols;
  sols = *eSols;

  if (debug) cout << "TtSemilepSolutionFilter: Found "<< sols.size()  << " Semilepton soutions.\n";

  //cout << "mass :";

  if (sols.size()!= 12)  return false;
  TtSemiEvtSolution & sol = sols[0];
 allSolution++;
  if ((sol.getCalLepb().pt()<jetPtCut) || (sol.getCalHadb().pt()<jetPtCut) ||  (sol.getCalHadp().pt()<jetPtCut) || (sol.getCalHadq().pt()<jetPtCut))
      return false;

  if ((sol.getCalLepb().et()<jetEtCut) || (sol.getCalHadb().et()<jetEtCut) ||  (sol.getCalHadp().et()<jetEtCut) || (sol.getCalHadq().et()<jetEtCut))
      return false;

  //FIXME: what is this in semilep case?
  //if (debug) cout << sol.getWpDecay() << " "<<sol.getWmDecay()<<endl;

  if (sol.getRecLepm().p4().pt()<leptonTriggerPtCut)
	return false;

  if (sol.getRecLepm().p4().pt()<leptonPtCut)
	return false;

  if (sol.getRecLepm().getTrackIso()>leptonTrackIsoCut)
	return false;
	
  if (sol.getRecLepm().getCaloIso()>leptonCaloIsoCut)
	return false;	

  //FIXME: do I need a similar cut? does this use MC info?
  // Check that we have an e-mu event

  //if ( ((sol.getWpDecay()!="muon")&&(sol.getWmDecay()!="muon")) ||
	//((sol.getWpDecay()!="electron")&& (sol.getWmDecay()!="electron")) )
	//return false;

//   if ( (sol.getWpDecay()=="electron")&&( !checkElectron(sol.getElectronp()) ) )
//      return false;
// 
//   if ( (sol.getWmDecay()=="electron")&&( !checkElectron(sol.getElectronm()) ) )
//      return false;

  ++goodSolution;
//FIXME:do I need this check?
   // try {if ( (abs(sols[0].getGenEvent()->lepton()->pdgId())==15) || 
  //(abs(sols[0].getGenEvent()->leptonBar()->pdgId())==15) ) ++tau;
   // } catch (...){exception++;}

//FIXME: can I use this?
   // if (sols[0].getJetB().getPartonFlavour()==5) ++B;
   //   else ++nonB;
   //if (sols[0].getJetBbar().getPartonFlavour()==5) ++B;
   //  else ++nonB;

  if (debug) cout << sol.getCalLepb().et() << " " <<sol.getCalHadb().et() << " " << sol.getCalHadp().et()<< " " << sol.getCalHadq().et()<< "" << endl; 
  if (debug) cout << sol.getCalLepb().pt() << " " <<sol.getCalHadb().pt() << " " << sol.getCalHadp().pt()<< " " << sol.getCalHadq().pt()<< "" << endl;
 if (debug) cout <<" ============================ End TtSemilepSolutionFilter ============================" <<endl;
  return true;
}


bool TtSemilepSolutionFilter::checkElectron (const TopElectron & electron) const
{
  if ( (electron.hadronicOverEm() > 0.1) ||
  	(electron.eSuperClusterOverP() > 3.) ||
	(electron.eSuperClusterOverP() < 0.8) ) return false;
  return true;

}

