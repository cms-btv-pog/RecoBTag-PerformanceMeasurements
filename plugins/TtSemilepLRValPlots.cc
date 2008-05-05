#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepLRValPlots.h"

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
TtSemilepLRValPlots::TtSemilepLRValPlots(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");
  rootFileName      = iConfig.getParameter<string> ("rootFileName");
  obsFileName      = iConfig.getParameter<string> ("obsFileName");  
  csa = iConfig.getParameter< bool > ("CSA");
  csaProcID = iConfig.getParameter< bool > ("CSAProcID");
  if (!csa) {
    weight = iConfig.getParameter< double > ("weight");
  }

  nrJetCombObs   = iConfig.getParameter<int> ("nrJetCombObs");
  obsNrs	   = iConfig.getParameter< vector<int> > ("JetCombObs");
  lrBins	   = iConfig.getParameter< int > ("lrBins");
  lrMin	   = iConfig.getParameter< double > ("lrMin");
  lrMax	   = iConfig.getParameter< double > ("lrMax");

  evtsols           = iConfig.getParameter<edm::InputTag> ("EvtSolution");
  lrFits = "pol4";

  bTagCutLabel = iConfig.getParameter< string >("bTagCutLabel");
  bCut = iConfig.getParameter< double >("bCut");

  theFile = new TFile("a.root", "RECREATE");
  theFile->cd();

  // define all histograms & fit functions
  myLRhelper = new LRHelpFunctions(lrBins, lrMin, lrMax, lrFits);
  myLRhelper->readObsHistsAndFits(obsFileName, obsNrs, false);

  goodSolutiontt0j = 0;
  goodSolutiontt1j = 0;
  goodSolutiontt2j = 0;
  goodSolutiontt3j = 0;
  goodSolutiontt4j = 0;
  goodSolutiontt0jother = 0;
  goodSolutiontt1jother = 0;
  goodSolutiontt2jother = 0;
  goodSolutiontt3jother = 0;
  goodSolutiontt4jother = 0;
  goodSolutionrest = 0;   
  goodSolutiontt0jbcut = 0;
  goodSolutiontt1jbcut = 0;
  goodSolutiontt2jbcut = 0;
  goodSolutiontt3jbcut = 0;
  goodSolutiontt4jbcut = 0;
  goodSolutiontt0jotherbcut = 0;
  goodSolutiontt1jotherbcut = 0;
  goodSolutiontt2jotherbcut = 0;
  goodSolutiontt3jotherbcut = 0;
  goodSolutiontt4jotherbcut = 0;
  goodSolutionrestbcut = 0; 
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
  
   cout << " Number of events before btag cut --------------------------------------+-+-+" << endl;  
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt0jbcut << "  tt0j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt1jbcut << "  tt1j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt2jbcut << "  tt2j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt3jbcut << "  tt3j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt4jbcut << "  tt4j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt0jotherbcut << "  tt0j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt1jotherbcut << "  tt1j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt2jotherbcut << "  tt2j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt3jotherbcut << "  tt3j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt4jotherbcut << "  tt4j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutionrestbcut << "  rest" << endl;
   cout << " Number of events after btag cut --------------------------------------+-+-+" << endl;  
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt0j << "  tt0j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt1j << "  tt1j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt2j << "  tt2j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt3j << "  tt3j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt4j << "  tt4j semimu" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt0jother << "  tt0j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt1jother << "  tt1j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt2jother << "  tt2j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt3jother << "  tt3j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutiontt4jother << "  tt4j other" << endl;
   cout << "TtSemiLepLRObsPlots: Events per 100pb-1 " << goodSolutionrest << "  rest" << endl;

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
    
    if (csa) { //reading out the event weights needed for processing the CSA07 soups
      edm::Handle< double> weightHandle;
      iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
      weight = * weightHandle;
      if (csaProcID) {
	int procID = csa07::csa07ProcessId(iEvent);
	if (debug) cout << "processID: " << procID
			<< " - name: " << csa07::csa07ProcessName(procID) 
			<< " - weight: "<< weight << endl;
      }
    }
    
    double bestlr = -999999.;
    int bestSol=-1;
    double maxProbChi2=-999;
    double maxbTagCut=-9999; 
    bool isSemiLeptonic=false;    
    int procID=0;

    //counting the number of events which past the cuts	
    for (int w = 0; w < 12; w++){  
      if (sols[w].getProbChi2()>maxProbChi2){maxProbChi2=sols[w].getProbChi2();}   
      if (sols[w].getCalHadb().getBDiscriminator(bTagCutLabel) > maxbTagCut) {maxbTagCut=sols[w].getCalHadb().getBDiscriminator(bTagCutLabel);} 
    }
   
    try {
      
      edm::Handle<TtGenEvent> genEvent;
      iEvent.getByLabel ("genEvt",genEvent);
      
      isSemiLeptonic = sols[0].getGenEvent()->isSemiLeptonic();
      //isSemiLeptonic = true;
      
    } catch (...){cout << "Exception\n";}

    //cout << "I stopped after the exception" << endl;

    //remark here should come some csa07 and geninfo booleans 
    if (csa) {
      if(maxProbChi2>0&&maxbTagCut>bCut) {    goodSolutionrest+=weight;}
      if (csaProcID){
	procID = csa07::csa07ProcessId(iEvent);
	if(debug) cout << maxProbChi2 << " " << maxbTagCut << endl;
	if(maxProbChi2>0&&maxbTagCut>bCut) {
	  if(procID==22){
	    if(isSemiLeptonic) {goodSolutiontt0j+=weight;} else {goodSolutiontt0jother+=weight;}
	  } else if (procID==23){
	    if(isSemiLeptonic) {goodSolutiontt1j+=weight;} else {goodSolutiontt1jother+=weight;}
	  } else if (procID==24){
	    if(isSemiLeptonic) {goodSolutiontt2j+=weight;} else {goodSolutiontt2jother+=weight;}
	  } else if (procID==25){
	    if(isSemiLeptonic) {goodSolutiontt3j+=weight;} else {goodSolutiontt3jother+=weight;}
	  } else if (procID==26){
	    if(isSemiLeptonic) {goodSolutiontt4j+=weight;} else {goodSolutiontt4jother+=weight;}
	  } else {goodSolutionrest+=weight;
	  }
	  if (debug) cout << "this event is selected" <<endl;
	}
	if(maxProbChi2>0) {
	  if(procID==22){
	    if(isSemiLeptonic) {goodSolutiontt0jbcut+=weight;} else {goodSolutiontt0jotherbcut+=weight;}
	  } else if (procID==23){
	    if(isSemiLeptonic) {goodSolutiontt1jbcut+=weight;} else {goodSolutiontt1jotherbcut+=weight;}
	  } else if (procID==24){
	    if(isSemiLeptonic) {goodSolutiontt2jbcut+=weight;} else {goodSolutiontt2jotherbcut+=weight;}
	  } else if (procID==25){
	    if(isSemiLeptonic) {goodSolutiontt3jbcut+=weight;} else {goodSolutiontt3jotherbcut+=weight;}
	  } else if (procID==26){
	    if(isSemiLeptonic) {goodSolutiontt4jbcut+=weight;} else {goodSolutiontt4jotherbcut+=weight;}
	  } else {goodSolutionrestbcut+=weight;
	  }
	}
      }
    }

    for (int i = 0; i < 12; ++i) {
      vector<double> obsVals;
      for(int j = 0; j < nrJetCombObs; j++){
	if( myLRhelper->obsFitIncluded(obsNrs[j]) )
          obsVals.push_back(sols[i].getLRJetCombObsVal(obsNrs[j]));
	//cout <<sols[i].getLRJetCombObsVal(obsNrs[j]) << endl;
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
    
    if (debug) cout << "Best LR value : " << bestlr << "  partonflavour : " << fabs(sols[bestSol].getCalLepb().getPartonFlavour()) <<endl;  
 
    try {

      edm::Handle<TtGenEvent> genEvent;
      iEvent.getByLabel ("genEvt",genEvent);
      
      matchB = (fabs(sols[bestSol].getCalLepb().getPartonFlavour())==5) ;	
      
    } catch (...){cout << "Exception\n";}
    
    if (debug) cout << "Best LR value" << bestlr <<endl;   
    if(maxProbChi2>0&&maxbTagCut>bCut){
      if (matchB) myLRhelper -> fillLRSignalHist(bestlr, weight);
      else myLRhelper -> fillLRBackgroundHist(bestlr, weight);
    }
  }
}
