#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepLRObsPlots.h"

#include "PhysicsTools/Utilities/interface/DeltaR.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtDilepLRObsPlots::TtDilepLRObsPlots(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");
  rootFileName      = iConfig.getParameter<string> ("rootFileName");
  weight = iConfig.getParameter< double > ("weight");

  nrSignalSelObs   = iConfig.getParameter<int> ("nrSignalSelObs");
  obsNrs	   = iConfig.getParameter< vector<int> > ("SignalSelObs");
  nrSignalSelHistBins	   = iConfig.getParameter< int > ("nrSignalSelHistBins");
  obsMin	   = iConfig.getParameter< vector<double> > ("SignalSelObsMin");
  obsMax	   = iConfig.getParameter< vector<double> > ("SignalSelObsMax");

  evtsols    = iConfig.getParameter<edm::InputTag>("evtSolution");
  jetSource_ = iConfig.getParameter<edm::InputTag>("jetSource");

  for(int j = 0; j < nrSignalSelObs; j++){
    obsFits.push_back("pol4");
  }
  
  theFile = new TFile("a.root", "RECREATE");
  theFile->cd();

  // define all histograms & fit functions
  myLRhelper = new LRHelpFunctions(obsNrs, nrSignalSelHistBins, obsMin, obsMax, obsFits);  

  allSolution = 0;
  goodSolution = 0;
  B=0;nonB=0;tau=0;

  myLRSignalSelObservables = new TtDilepLRSignalSelObservables();
  myLRSignalSelObservables->jetSource(jetSource_);
}


TtDilepLRObsPlots::~TtDilepLRObsPlots()
{
   delete myLRhelper;
   theFile->Close();
}

void TtDilepLRObsPlots::endJob()
{
   cout << "TtDilepLRObsPlots: Events " << goodSolution <<endl;
   cout << "TtDilepLRObsPlots: Events per fb-1 " << goodSolution*weight << endl;
   pair<double, double> match = myLRhelper->getBnumbers();
cout << "TtDilepLRObsPlots match " << B*weight<<" "<<nonB*weight<< " "<< tau<<" "<<match.second<<endl;

   // store histograms and fits in root-file
   myLRhelper -> storeToROOTfile(rootFileName);
   if (debug) cout << "************* Finished writing histograms to file" << endl;

}

void
TtDilepLRObsPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   // get the event solution
  edm::Handle< std::vector<TtDilepEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);

  std::vector<TtDilepEvtSolution> sols;
  sols = *eSols;
  if (debug) cout << "TtDilepLRObsPlots: Found "<< sols.size()  << " dilepton soutions.\n";

  if(sols.size()== 2) {

    ++goodSolution;

    vector < vector< TtDilepLRSignalSelObservables::IntBoolPair > > obsMatch;
    vector <int> matchSum;
    int bestSol=-1, bestMatch = -1;

    for (int i = 0; i < 2; ++i) {
      obsMatch.push_back(myLRSignalSelObservables->operator()(sols[i],iEvent, true));
    }
    for (int i = 0; i < 2; ++i) {
      // Fill the observables 
      for(int j = 0; j < nrSignalSelObs; ++j) {
	if( myLRhelper->obsFitIncluded(obsNrs[j]) ) {
          if (obsMatch[i][obsNrs[j]-1].second) {
	    myLRhelper->fillToSignalHists(obsNrs[j],
		  sols[i].getLRSignalEvtObsVal(obsNrs[j]),
		  weight/(obsMatch[0][obsNrs[j]-1].second+obsMatch[1][obsNrs[j]-1].second));
	    for(int k = j; k < nrSignalSelObs; ++k) {
              if (obsMatch[i][obsNrs[k]-1].second) {
		myLRhelper->fillToSignalCorrelation(
		  obsNrs[j], sols[i].getLRSignalEvtObsVal(obsNrs[j]),
		  obsNrs[k], sols[i].getLRSignalEvtObsVal(obsNrs[k]),
		  weight/(obsMatch[0][obsNrs[j]-1].second*obsMatch[0][obsNrs[k]-1].second+
		  	obsMatch[1][obsNrs[j]-1].second*obsMatch[1][obsNrs[k]-1].second));
// cout << "Fill :" << obsMatch[0][obsNrs[j]-1].second*obsMatch[0][obsNrs[k]-1].second<<" "<<
// 		  	obsMatch[0][obsNrs[j]-1].second<<" "<<obsMatch[0][obsNrs[k]-1].second<<" "<<
// 			obsNrs[j]<<" "<< sols[i].getLRSignalEvtObsVal(obsNrs[j])<<" "<<
// obsNrs[k]<<" "<< sols[i].getLRSignalEvtObsVal(obsNrs[k])<<endl;
	      }
	    }
	  } else myLRhelper->fillToBackgroundHists(obsNrs[j],
		  sols[i].getLRSignalEvtObsVal(obsNrs[j]), weight/(!obsMatch[0][obsNrs[j]-1].second+!obsMatch[1][obsNrs[j]-1].second));
	}
	
      }


      int matchSum = 0;
      for(int j = 0; j < nrSignalSelObs; j++){
	if (obsMatch[i][j].second) ++matchSum;
	if (debug) cout <<obsMatch[i][j].second;
      }
      if (debug) cout <<endl;
      if (matchSum>bestMatch) {
        bestMatch = matchSum;
	bestSol = i;
      }
    }
    if (debug) cout << "\nBest Solution:" <<sols[0].getBestSol() << " - "<<sols[1].getBestSol()<< " = "<<bestMatch<<" = "<<bestSol<< endl;

    B+=obsMatch[bestSol][0].second;
    B+=obsMatch[bestSol][1].second;
    nonB+=(!obsMatch[bestSol][0].second);
    nonB+=(!obsMatch[bestSol][1].second);
//  cout<< B<<nonB<< obsMatch[bestSol][0].second<<obsMatch[bestSol][1].second<<endl;
    if (sols[bestSol].getJetB().getPartonFlavour()==5) ++tau;
    if (sols[bestSol].getJetBbar().getPartonFlavour()==5) ++tau;

  }
  if (debug) cout <<" ============================ End TtDilepLRObsPlots ============================" <<endl;
}


bool TtDilepLRObsPlots::allMatch(const TtDilepEvtSolution & sol) 
{
  if (debug) cout<<"match test " << endl;
  if (debug) cout << sol.getJetB().getPartonFlavour()<< " "<<sol.getJetBbar().getPartonFlavour()<< " ";
  if (debug) cout << sol.getWpDecay() << " "<<sol.getWmDecay()<<endl;
// cout << sol.getCalJetB().et()<< " " << sol.getCalJetBbar().et()<<endl;
// cout << sol.getGenB()->et()<< " " << sol.getGenBbar()->et()<<endl;

  if (debug) cout << "rec pT/eta/phi: "<< sol.getJetB().pt()<< " / " << sol.getJetB().eta()<< " / " <<sol.getJetB().phi()
    << " -- " <<sol.getJetBbar().pt()<< " / " << sol.getJetBbar().eta()<< " / " <<sol.getJetBbar().phi()
    << endl;
  if (debug) cout << "deltaR: "<< DeltaR<reco::Particle>()(sol.getCalJetB(), *(sol.getGenB()))
    << " " <<DeltaR<reco::Particle>()(sol.getCalJetBbar(), *(sol.getGenBbar()))
    << endl;
  if (debug) cout << "gen pT/eta/phi: "<< sol.getGenB()->pt()<< " / " << sol.getGenB()->eta()<< " / " <<sol.getGenB()->phi()
    << " -- " <<sol.getGenBbar()->pt()<< " / " << sol.getGenBbar()->eta()<< " / " <<sol.getGenBbar()->phi()
    << endl;
//   if (DeltaR<reco::Particle>()(sol.getCalJetB(), *(sol.getGenB()))>0.5) return false;
//   if (DeltaR<reco::Particle>()(sol.getCalJetBbar(), *(sol.getGenBbar()))>0.5) return false;

 if ((sol.getJetB().getPartonFlavour()!=5) || 
 (sol.getJetBbar().getPartonFlavour()!=5)) return false;

  if (debug) cout <<"getGenLepp() "<<	 sol.getGenEvent()->leptonBar()->pdgId()
		  <<" - getGenLepm() "<< sol.getGenEvent()->lepton()->pdgId()<<endl;

  bool lept = false;
  if ( ((sol.getWpDecay()=="muon")&&(abs(sol.getGenEvent()->leptonBar()->pdgId())==13)) &&
   ((sol.getWmDecay()=="electron")&&(abs(sol.getGenEvent()->lepton()->pdgId())==11)) ) lept = true;

  if ( ((sol.getWpDecay()=="electron")&&(abs(sol.getGenEvent()->leptonBar()->pdgId())==11)) &&
   ((sol.getWmDecay()=="muon")&&(abs(sol.getGenEvent()->lepton()->pdgId())==13)) ) lept = true;
  if ( (abs(sol.getGenEvent()->lepton()->pdgId())==15) || 
  (abs(sol.getGenEvent()->leptonBar()->pdgId())==15) ) ++tau;
  if (debug) cout << "Lepton: "<< lept << " " 
    << DeltaR<reco::Particle>()(sol.getLeptNeg(), *(sol.getGenLepm()))
    << " " << DeltaR<reco::Particle>()(sol.getLeptPos(), *(sol.getGenLepp()))
    <<endl;
  if (!lept) return false;
  if (DeltaR<reco::Particle>()(sol.getLeptNeg(), *(sol.getGenLepm())) > 0.1) return false;
  if (DeltaR<reco::Particle>()(sol.getLeptPos(), *(sol.getGenLepp())) > 0.1) return false;
  if (debug) cout << "Match OK\n";
  return true;

}
