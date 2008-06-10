#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepBTagAnalysis.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtSemilepBTagAnalysis::TtSemilepBTagAnalysis(const edm::ParameterSet& iConfig)
{
  if (debug) cout << "TtSemilepBTagAnalysis ctor" << endl;
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");

  rootFileName      = iConfig.getParameter<string> ("rootFileName");

  obsFileName      = iConfig.getParameter<string> ("obsFileName");
  nrJetCombObs   = iConfig.getParameter<int> ("nrJetCombObs");
  obsNrs	   = iConfig.getParameter< vector<int> > ("JetCombObs");

  evtsols           = iConfig.getParameter<edm::InputTag> ("EvtSolution");
  weight = iConfig.getParameter< double > ("weight");



  bool update = iConfig.getParameter<bool>( "update" );

  vector<edm::ParameterSet> taggerConfig = iConfig.getParameter< vector<edm::ParameterSet> >("tagConfig");

  // eta and pt ranges (bins for differential plots)
  etaRanges = iConfig.getParameter< vector<double> >("etaRanges");
  etRanges = iConfig.getParameter< vector<double> >("etRanges");

  lrBins	   = iConfig.getParameter< int > ("lrBins");
  lrMin	   = iConfig.getParameter< double > ("lrMin");
  lrMax	   = iConfig.getParameter< double > ("lrMax");

  mcMode = iConfig.getParameter<bool>( "mcMode" );


  if (update) {
    TString inputFileName = iConfig.getParameter<string> ("inputFileName");
    theFile = new TFile (inputFileName) ;
    if (debug) cout << "Open file in upate mode: " << inputFileName<<endl;
  } else {
    theFile = new TFile (TString (rootFileName) , "RECREATE" ) ;
    theFile->cd();
  }
  TH1::AddDirectory(kFALSE);

  if (debug) cout << "LHR plots" << endl;
  // define all histograms & fit functions
  //myLRhelper = new LRHelpFunctions(lrBins, lrMin, lrMax, lrFits);
  myLRhelper = new LRHelpFunctions();
  if (debug) cout << "Read LHR observables" << endl;
  myLRhelper->readObsHistsAndFits(obsFileName, obsNrs, false);
  if (debug) cout << "LHR observables read" << endl;

  if (debug) cout << "Create analysis histos" << endl;
  for (unsigned int iWP = 0; iWP != taggerConfig.size(); ++iWP) {

    string bTaggerName = taggerConfig[iWP].getParameter<string>("tagger");
    double cut = taggerConfig[iWP].getParameter<double>("cut");
    bTagger.push_back( StringIntPair(bTaggerName, cut) );

    TString baseName(bTaggerName);
    baseName+="_";
    baseName+=cut;
    baseName.ReplaceAll ( " " , "" );

    if (debug) cout << "Book plots for Working Point " << baseName<<endl;

    taggedJetsHistos.push_back( new TtEtEtaHistoCollector("TaggedJets_"+baseName, etRanges, etaRanges,
	lrMin, lrMax, lrBins, update) ) ;
    taggedJetsHistos.back()->Sumw2();
    if (mcMode) {
      taggedLightJetsHistos.push_back( new TtEtEtaHistoCollector("TaggedLightJets_"+baseName, etRanges, etaRanges,
	lrMin, lrMax, lrBins, update) ) ;
      taggedLightJetsHistos.back()->Sumw2();
      taggedBJetsHistos.push_back( new TtEtEtaHistoCollector("TaggedBJets_"+baseName, etRanges, etaRanges,
	lrMin, lrMax, lrBins, update) ) ;
      taggedBJetsHistos.back()->Sumw2();
    }

  }

  if (debug) cout << "Book plots for AllEvents" << endl;

  TString baseName = TString("AllJets");
  allJetsHistos = new TtEtEtaHistoCollector(baseName, etRanges, etaRanges,
      lrMin, lrMax, lrBins, update);
  allJetsHistos->Sumw2();

  if (mcMode) {
    if (debug) cout << "Monte Carlo mode\n";
    baseName = TString("LRsignal");
    lrSHistos = new TtEtEtaHistoCollector(baseName, etRanges, etaRanges,
	lrMin, lrMax, lrBins, update);
    lrSHistos->Sumw2();

    baseName = TString("LRbackground");
    lrBHistos = new TtEtEtaHistoCollector(baseName, etRanges, etaRanges,
	lrMin, lrMax, lrBins, update);
    lrBHistos->Sumw2();

    baseName = TString("All_LightJets");
    allLightHistos = new TtEtEtaHistoCollector(baseName, etRanges, etaRanges,
	lrMin, lrMax, lrBins, update);
    allLightHistos->Sumw2();

  }
  
  if (update) {
    theFile->Close();
  }

  if (debug) cout << "All Done! " << endl;

//   // eta jet
//   double etaMin = iConfig.getParameter<double>("etaMin");
//   double etaMax = iConfig.getParameter<double>("etaMax");
//   // et jet
//   double etMin = iConfig.getParameter<double>("etMin");
//   double etMax = iConfig.getParameter<double>("etMax");
// 
//   // specify jet and parton kinematic cuts.
//   jetSelector.setEtaMin(etaMin);
//   jetSelector.setEtaMax(etaMax);
//   jetSelector.setEtMin(etaMin);
//   jetSelector.setEtMax(etaMax);
}


TtSemilepBTagAnalysis::~TtSemilepBTagAnalysis()
{
  if (debug) cout << "TtSemilepBTagAnalysis Destructor" << endl;
  delete myLRhelper;
  theFile->Close();
  delete allJetsHistos;
  if (mcMode) {
    delete lrSHistos;
    delete lrBHistos;
    delete allLightHistos;
  }

  for (unsigned int iWP = 0; iWP != bTagger.size(); ++iWP) {
    delete taggedJetsHistos[iWP];
    if (mcMode) delete taggedLightJetsHistos[iWP];
    if (mcMode) delete taggedBJetsHistos[iWP];
  }
}

void TtSemilepBTagAnalysis::writePlotter(TtEtEtaHistoCollector* plotter)
{
//  theFile->mkdir(plotter->baseName());
//  theFile->cd(plotter->baseName());
  plotter->write();
//  gFile->cd();
}

void TtSemilepBTagAnalysis::endJob()
{
    theFile = new TFile (TString (rootFileName) , "RECREATE" ) ;
    theFile->cd();

  if (debug) cout << "TtSemilepBTagAnalysis endJob" << endl;
  writePlotter(allJetsHistos);
  if (mcMode) {
    writePlotter(lrSHistos);
    writePlotter(lrBHistos);
    writePlotter(allLightHistos);
  }
  for (unsigned int iWP = 0; iWP != bTagger.size(); ++iWP) {
    writePlotter(taggedJetsHistos[iWP]);
    if (mcMode) writePlotter(taggedLightJetsHistos[iWP]);
    if (mcMode) writePlotter(taggedBJetsHistos[iWP]);
  }

//   TtEtEtaHistoCollector allIntJetsHistos = allJetsHistos->buildIntegratedCollector();
//   writePlotter(&allIntJetsHistos);
// 
//   if (mcMode) {
//     TtEtEtaHistoCollector intLrSHistos  =  lrSHistos->buildIntegratedCollector();
//     TtEtEtaHistoCollector intLrBHistos  =  lrBHistos->buildIntegratedCollector();
//     TtEtEtaHistoCollector intallLightHistos = allLightHistos->buildIntegratedCollector();
//     writePlotter(&intLrSHistos );
//     writePlotter(&intLrBHistos );
//     writePlotter(&intallLightHistos);
// 
//     TtEtEtaHistoCollector xb = intLrSHistos.buildRatioCollector(TString("B_Purity"),
// 			allIntJetsHistos);
//     writePlotter(&xb);
// 
//   }
  
//   for (unsigned int iWP = 0; iWP != bTagger.size(); ++iWP) {
//   TtEtEtaHistoCollector intLrSHistos  =  taggedJetsHistos[iWP]->buildIntegratedCollector();
//     writePlotter(taggedJetsHistos[iWP]);
//     if (mcMode) writePlotter(taggedLightJetsHistos[iWP]);
//   }


}

void
TtSemilepBTagAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (debug) cout << "************* TtSemilepBTagAnalysis Start" << endl;
  // get the event solution
  edm::Handle< std::vector<TtSemiEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);

  std::vector<TtSemiEvtSolution> sols;
  sols = *eSols;
 
  //reading out the event weights needed for processing the CSA07 soups
    
  edm::Handle< double> weightHandle;
  iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
  double CSA07weight = * weightHandle; 
  //this line will overwrite the weight as stated in the .cfi file
  weight = CSA07weight;
  
  if(sols.size()== 12) {
    
    //vector < double > lr;
    double bestlr = -999999.;
    int bestSol=-1;
    //double highestLR = -10000000000000.;
    double maxProbChi2=-999;
    double maxbTagCut=-9999; 
    
    //defining extra cuts	
    for (int w = 0; w < 12; w++){  
      if (sols[w].getProbChi2()>maxProbChi2){maxProbChi2=sols[w].getProbChi2();}   
      if (sols[w].getCalHadb().getBDiscriminator(bTagCutLabel) > maxbTagCut) {maxbTagCut=sols[w].getCalHadb().getBDiscriminator(bTagCutLabel);} 
    }
    
    //calculate the combined lr value, save the highest one and the corresponding solution
    for (int i = 0; i < 12; ++i) {
      vector<double> obsVals;
      for(int j = 0; j < nrJetCombObs; j++){
	if( myLRhelper->obsFitIncluded(obsNrs[j]) ) obsVals.push_back(sols[i].getLRJetCombObsVal(obsNrs[j]));
	// cout << "Obs " << j<<" "<<sols[i].getLRSignalEvtObsVal(obsNrs[j])<<endl;
      }
      if (debug) cout << "start calculating Combined LR value" << endl;
      double lr = myLRhelper->calcLRval(obsVals);
      if (debug) cout << "Combined LR value : " << lr << endl;
      if (lr >= bestlr) {
	bestlr = lr;
	bestSol = i;
      }
    } 
    
    bool matchB = false;
    
    if (debug) cout << "Best LR value : " << bestlr << "  partonflavour : " << fabs(sols[bestSol].getCalLepb().getPartonFlavour()) <<endl;  
 
    try {

      edm::Handle<TtGenEvent> genEvent;
      iEvent.getByLabel ("genEvt",genEvent);
      
      matchB = (fabs(sols[bestSol].getCalLepb().getPartonFlavour())==5) ;	
      
    } catch (...){cout << "Exception\n";}       
    
    if(maxProbChi2>0&&maxbTagCut>bCut){
      analyzeJet(sols[bestSol].getCalLepb(), bestlr, matchB);
    }
    // analyzeJet(sols[bestSol].getCalJetBbar(), bestlr, matchBbar);
  }
}

void TtSemilepBTagAnalysis::analyzeJet(const TopJet& jet, double lr, bool matchB)
{
//   int flavour = jet.getPartonFlavour();
//   matchB = flavour==5;

  allJetsHistos->analyze(jet, lr, weight);

  if (mcMode) {
    if (matchB) lrSHistos->analyze(jet, lr, weight);
      else lrBHistos->analyze(jet, lr, weight);
    if (!matchB) {
      allLightHistos->analyze(jet, lr, weight);
      if (debug) cout << "Light: " <<jet.getPartonFlavour()<<" ";
    }
  }

  for (unsigned int iWP = 0; iWP != bTagger.size(); ++iWP) {
    if (debug) cout << "B-Jet: " << bTagger[iWP].first<< " -> " << jet.getBDiscriminator(bTagger[iWP].first) <<endl;
    if (jet.getBDiscriminator(bTagger[iWP].first) > bTagger[iWP].second ) {
      taggedJetsHistos[iWP]->analyze(jet, lr, weight);
      if (mcMode) {
	if (!matchB) taggedLightJetsHistos[iWP]->analyze(jet, lr, weight);
	if (matchB) taggedBJetsHistos[iWP]->analyze(jet, lr, weight);
      }
    }
    if (debug) if (!matchB) {
      if (jet.getBDiscriminator(bTagger[iWP].first) > bTagger[iWP].second ) cout <<"tagged\n";
      else cout <<"NOT\n";
    }
  }

}
