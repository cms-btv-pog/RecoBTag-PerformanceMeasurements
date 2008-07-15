#include "RecoBTag/PerformanceMeasurements/plugins/TtDilepLRValPlots.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "CSA07EffAnalyser/CSA07EffAnalyser/interface/CSA07ProcessId.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtDilepLRValPlots::TtDilepLRValPlots(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");
  rootFileName      = iConfig.getParameter<string> ("rootFileName");
  obsFileName      = iConfig.getParameter<string> ("obsFileName");
  csa = iConfig.getParameter< bool > ("CSA");
  if (!csa) {
    weight = iConfig.getParameter< double > ("weight");
  }

  nrSignalSelObs   = iConfig.getParameter<int> ("nrSignalSelObs");
  obsNrs	   = iConfig.getParameter< vector<int> > ("SignalSelObs");
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


TtDilepLRValPlots::~TtDilepLRValPlots()
{
   if (debug) cout << "TtDilepLRValPlots Destructor" << endl;
   delete myLRhelper;
   theFile->Close();



}
void TtDilepLRValPlots::endJob()
{
   if (debug) cout << "TtDilepLRValPlots endJob" << endl;
   // store histograms and fits in root-file
   myLRhelper->makeAndFitPurityHists();
   myLRhelper->storeToROOTfile(rootFileName);
}

void
TtDilepLRValPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (debug) cout << "************* TtDilepLRValPlots Start" << endl;
  // get the event solution
  edm::Handle< std::vector<TtDilepEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);

  std::vector<TtDilepEvtSolution> sols;
  sols = *eSols;

  if(sols.size()== 2) {

    if (csa) {
      edm::Handle< double> weightHandle;
      iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
      weight = * weightHandle;
      int procID = csa07::csa07ProcessId(iEvent);
      if (debug) cout << "processID: " << procID
	<< " - name: " << csa07::csa07ProcessName(procID) 
	<< " - weight: "<< weight << endl;
    }

    vector < double > lr;
    int bestSol=-1;
    //double highestLR = -10000000000000.;

    for (int i = 0; i < 2; ++i) {
      vector<double> obsVals;
      for(int j = 0; j < nrSignalSelObs; j++){
	if( myLRhelper->obsFitIncluded(obsNrs[j]) )
          obsVals.push_back(sols[i].getLRSignalEvtObsVal(obsNrs[j]));
//	  cout << "Obs " << j<<" "<<sols[i].getLRSignalEvtObsVal(obsNrs[j])<<endl;
//       cout << j<<myLRhelper->obsFitIncluded(obsNrs[j])<<" ";
      }
//       cout << endl;
      // Fill the LR 
      lr.push_back(myLRhelper->calcPtdrLRval(obsVals) );

//       if (lr[i]>highestLR) {
//         highestLR = lr[i];
// 	bestSol = i;
//       }
    }
    if (lr[0]>lr[1]) bestSol=0;
      else bestSol=1;
//  cout << "TtDilepBTagAnalysis: " << lr[0]<<" "<<lr[1]<<" "<<bestSol<<endl;

    //if (debug)
//     if (lr[bestSol] < -10) 
//     cout << "Best Solution:" <<lr[0] << " - "<<lr[1]<< " = "
//     <<bestSol<<" = "<<sols[bestSol].getJetB().getPartonFlavour()<<
//     sols[bestSol].getJetBbar().getPartonFlavour()<< endl;

  bool matchB = false;
  bool matchBbar = false;

  try {
    // cout <<endl;
    double dr1, dr2;
    edm::Handle<TtGenEvent> genEvent;
    iEvent.getByLabel ("genEvt",genEvent);

    if (genEvent->numberOfBQuarks()>1) {
      dr1 = DeltaR<reco::Particle>()(sols[bestSol].getCalJetB(), *(genEvent->b()));
      dr2 = DeltaR<reco::Particle>()(sols[bestSol].getCalJetB(), *(genEvent->bBar()));
//       matchB = ( (solution.getJetB().getPartonFlavour()==5) && (dr<0.4) );
//       matchB = ( (dr1<0.4) || (dr2<0.4));
      matchB = ( (dr1<0.4) || (dr2<0.4) || (sols[bestSol].getCalJetB().getPartonFlavour()==5));

//       cout <<"test "<< matchB <<" "<< DeltaR<reco::Particle>()(sols[bestSol].getCalJetB(), *(genEvent->b())) 
//       <<" "<< DeltaR<reco::Particle>()(sols[0].getCalJetB(), *(genEvent->b())) 
//     <<" "<< DeltaR<reco::Particle>()(sols[1].getCalJetB(), *(genEvent->b())) <<endl;

      dr1 = DeltaR<reco::Particle>()(sols[bestSol].getCalJetBbar(), *(genEvent->b()));
      dr2 = DeltaR<reco::Particle>()(sols[bestSol].getCalJetBbar(), *(genEvent->bBar()));
      matchBbar = ( ( (dr1<0.4) || (dr2<0.4) || (sols[bestSol].getCalJetBbar().getPartonFlavour()==5)) );
//        matchBbar = ( (solution.getJetBbar().getPartonFlavour()==5) && (dr<0.4) );
//       matchBbar = ( ( (dr1<0.4) || (dr2<0.4)) );
    }
  } catch (...){cout << "Exception\n";}

    if (matchB) myLRhelper -> fillLRSignalHist(lr[bestSol], weight);
      else myLRhelper -> fillLRBackgroundHist(lr[bestSol], weight);
    if (matchBbar) myLRhelper -> fillLRSignalHist(lr[bestSol], weight);
      else myLRhelper -> fillLRBackgroundHist(lr[bestSol], weight);
  }
}
