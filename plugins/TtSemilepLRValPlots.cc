#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepLRValPlots.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtSemilepLRValPlots::TtSemilepLRValPlots(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");
  rootFileName      = iConfig.getParameter<string> ("rootFileName");
  obsFileName      = iConfig.getParameter<string> ("obsFileName");
  weight = iConfig.getParameter< double > ("weight");

  nrJetCombObs   = iConfig.getParameter<int> ("nrJetCombObs");
  obsNrs	   = iConfig.getParameter< vector<int> > ("JetCombObs");
  lrBins	   = iConfig.getParameter< int > ("lrBins");
  lrMin	   = iConfig.getParameter< double > ("lrMin");
  lrMax	   = iConfig.getParameter< double > ("lrMax");

  evtsols           = iConfig.getParameter<edm::InputTag> ("EvtSolution");
  lrFits = "pol4";

  theFile = new TFile("a.root", "RECREATE");
  theFile->cd();

  // define all histograms & fit functions
  myLRhelper = new LRHelpFunctions(lrBins, lrMin, lrMax, lrFits);
  myLRhelper->readObsHistsAndFits(obsFileName, obsNrs, false);
}


TtSemilepLRValPlots::~TtSemilepLRValPlots()
{
   if (debug) cout << "TtSemilepLRValPlots Destructor" << endl;
   delete myLRhelper;
   theFile->Close();



}
void TtSemilepLRValPlots::endJob()
{
   if (debug) cout << "TtSemilepLRValPlots endJob" << endl;
   // store histograms and fits in root-file
   myLRhelper->makeAndFitPurityHists();
   myLRhelper->storeToROOTfile(rootFileName);
}

void
TtSemilepLRValPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (debug) cout << "************* TtSemilepLRValPlots Start" << endl;
  // get the event solution
  edm::Handle< std::vector<TtSemiEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);
  
  std::vector<TtSemiEvtSolution> sols;
  sols = *eSols;
  
  if(sols.size()== 12) {
    
    //vector < double > lr;
    double bestlr = -999999.;
    int bestSol=-1;
    //double highestLR = -10000000000000.;
    
    for (int i = 0; i < 12; ++i) {
      vector<double> obsVals;
      for(int j = 0; j < nrJetCombObs; j++){
	if( myLRhelper->obsFitIncluded(obsNrs[j]) )
          obsVals.push_back(sols[i].getLRJetCombObsVal(obsNrs[j]));
	cout <<sols[i].getLRJetCombObsVal(obsNrs[j]) << endl;
	//	  cout << "Obs " << j<<" "<<sols[i].getLRSignalEvtObsVal(obsNrs[j])<<endl;
	//       cout << j<<myLRhelper->obsFitIncluded(obsNrs[j])<<" ";
      }
      //       cout << endl;
      // Fill the LR 
      //FIXME: this should not be the one from ptdr
      //lr.push_back(myLRhelper->calcPtdrLRval(obsVals) );
      //double lr = myLRhelper->calcPtdrLRval(obsVals);
      if (debug) cout << "start calculating Combined LR value" << endl;
      double lr = myLRhelper->calcLRval(obsVals);
      if (debug) cout << "Combined LR value : " << lr << endl;
      if (lr >= bestlr) {
	bestlr = lr;
	bestSol = i;
      }
      //       if (lr[i]>highestLR) {
      //         highestLR = lr[i];
      // 	bestSol = i;
      //       }
    } 

    //else bestSol=0;
    //  cout << "TtSemilepBTagAnalysis: " << lr[0]<<" "<<lr[1]<<" "<<bestSol<<endl;
    
    //if (debug)
    //     if (lr[bestSol] < -10) 
    //     cout << "Best Solution:" <<lr[0] << " - "<<lr[1]<< " = "
    //     <<bestSol<<" = "<<sols[bestSol].getJetB().getPartonFlavour()<<
    //     sols[bestSol].getJetBbar().getPartonFlavour()<< endl;
    
    bool matchB = false;
     
    try {

      edm::Handle<TtGenEvent> genEvent;
      iEvent.getByLabel ("genEvt",genEvent);
      
      matchB = (sols[bestSol].getCalLepb().getPartonFlavour()==5) ;	
      
    } catch (...){cout << "Exception\n";}
    
    if (debug) cout << "Best LR value" << bestlr <<endl;
    if (matchB) myLRhelper -> fillLRSignalHist(bestlr, weight);
    else myLRhelper -> fillLRBackgroundHist(bestlr, weight);
   
  }
}
