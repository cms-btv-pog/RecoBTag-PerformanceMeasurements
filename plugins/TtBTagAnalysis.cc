#include "RecoBTag/PerformanceMeasurements/plugins/TtBTagAnalysis.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtBTagAnalysis::TtBTagAnalysis(const edm::ParameterSet& iConfig)
{
  if (debug) cout << "TtBTagAnalysis ctor" << endl;
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");

  rootFileName      = iConfig.getParameter<string> ("rootFileName");

  obsFileName      = iConfig.getParameter<string> ("obsFileName");
  nrSignalSelObs   = iConfig.getParameter<int> ("nrSignalSelObs");
  obsNrs	   = iConfig.getParameter< vector<int> > ("SignalSelObs");

  evtsols           = iConfig.getParameter<edm::InputTag> ("EvtSolution");
  weight = iConfig.getParameter< double > ("weight");

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



  if (debug) cout << "LHR plots" << endl;
  // define all histograms & fit functions
  myLRhelper = new LRHelpFunctions();
  if (debug) cout << "Read LHR observables" << endl;
  myLRhelper->readObsHistsAndFits(obsFileName, obsNrs, false);
  if (debug) cout << "LHR observables read" << endl;

  if (update) {
    TString inputFileName = iConfig.getParameter<string> ("inputFileName");
    theFile = new TFile (inputFileName) ;
    if (debug) cout << "Open file in upate mode: " << inputFileName<<endl;
  } else {
    theFile = new TFile (TString (rootFileName) , "RECREATE" ) ;
    theFile->cd();
    if (debug) cout << "Open file: " << rootFileName<<endl;
  }
  TH1::AddDirectory(kFALSE);
  if (debug) cout << "Create analysis histos" << endl;
  for (unsigned int iWP = 0; iWP != taggerConfig.size(); ++iWP) {

    string bTaggerName = taggerConfig[iWP].getParameter<string>("tagger");
    double cut = taggerConfig[iWP].getParameter<double>("cut");
    string baseName = taggerConfig[iWP].getParameter<string>("histoName");
    bTagger.push_back( StringIntPair(bTaggerName, cut) );

//     TString baseName(bTaggerName);
//     baseName+="_";
//     baseName+=cut;
//     baseName.ReplaceAll ( " " , "" );

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


TtBTagAnalysis::~TtBTagAnalysis()
{
  if (debug) cout << "TtBTagAnalysis Destructor" << endl;
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

void TtBTagAnalysis::writePlotter(TtEtEtaHistoCollector* plotter)
{
//  theFile->mkdir(plotter->baseName());
//  theFile->cd(plotter->baseName());
  plotter->write();
//  gFile->cd();
}

void TtBTagAnalysis::endJob()
{
    theFile = new TFile (TString (rootFileName) , "RECREATE" ) ;
    theFile->cd();

  if (debug) cout << "TtBTagAnalysis endJob" << endl;
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
TtBTagAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (debug) cout << "************* TtBTagAnalysis Start" << endl;
  // get the event solution
  edm::Handle< std::vector<TtDilepEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);

  std::vector<TtDilepEvtSolution> sols;
  sols = *eSols;

  if(sols.size()== 2) {

    vector < double > lr;
    int bestSol=-1;

    for (int i = 0; i < 2; ++i) {
      vector<double> obsVals;
      for(int j = 0; j < nrSignalSelObs; j++){
	if( myLRhelper->obsFitIncluded(obsNrs[j]) )
          obsVals.push_back(sols[i].getLRSignalEvtObsVal(obsNrs[j]));
//	  cout << "Obs " << j<<" "<<sols[i].getLRSignalEvtObsVal(obsNrs[j])<<endl;
      }
      lr.push_back(myLRhelper->calcPtdrLRval(obsVals) );
//       if (lr[i]>highestLR) {
//         highestLR = lr[i];
// 	bestSol = i;
//       }
    }
    if (lr[0]>lr[1]) bestSol=0;
      else bestSol=1;
//  cout << "TtDilepLRValPlots: " << lr[0]<<" "<<lr[1]<<" "<<bestSol<<endl;
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
//       matchB = ( (dr1<0.4) || (dr2<0.4) );
//       matchB = ( (sols[bestSol].getCalJetB().partonFlavour()==5));
      matchB = ( (dr1<0.4) || (dr2<0.4) || (sols[bestSol].getCalJetB().partonFlavour()==5));

      if (debug)  if ((!matchB)|| (sols[bestSol].getCalJetB().partonFlavour()!=5))
      cout <<"test "<< matchB<<" "<<sols[bestSol].getCalJetB().partonFlavour() <<" "<< DeltaR<reco::Particle>()(sols[bestSol].getCalJetB(), *(genEvent->b())) 
      <<" "<< dr1 <<" "<< dr2 <<endl;

      dr1 = DeltaR<reco::Particle>()(sols[bestSol].getCalJetBbar(), *(genEvent->b()));
      dr2 = DeltaR<reco::Particle>()(sols[bestSol].getCalJetBbar(), *(genEvent->bBar()));
      matchBbar = ( ( (dr1<0.4) || (dr2<0.4) || (sols[bestSol].getCalJetBbar().partonFlavour()==5)) );
//       matchBbar = ( (dr1<0.4) || (dr2<0.4) );
//       matchBbar = sols[bestSol].getCalJetBbar().partonFlavour()==5;

      if (debug) if ((!matchBbar)|| (sols[bestSol].getCalJetBbar().partonFlavour()!=5))
      cout <<"test "<< matchBbar<<" "<< sols[bestSol].getCalJetBbar().partonFlavour() <<" "<< DeltaR<reco::Particle>()(sols[bestSol].getCalJetBbar(), *(genEvent->bBar())) 
      <<" "<< dr1 <<" "<< dr2 <<endl;
    }
  } catch (...){cout << "Exception\n";}



    analyzeJet(sols[bestSol].getCalJetB(), lr[bestSol], weight, matchB);
    analyzeJet(sols[bestSol].getCalJetBbar(), lr[bestSol], weight, matchBbar);
  }
}

void TtBTagAnalysis::analyzeJet(const pat::Jet& jet, double lr, double weight, bool matchB)
{
//   int flavour = jet.partonFlavour`();
//   matchB = flavour==5;

  allJetsHistos->analyze(jet, lr, weight);

  if (mcMode) {
    if (matchB) lrSHistos->analyze(jet, lr, weight);
      else lrBHistos->analyze(jet, lr, weight);
    if (!matchB) {
      allLightHistos->analyze(jet, lr, weight);
      if (debug) cout << "Light: " <<jet.partonFlavour()<<" ";
    }
  }

  for (unsigned int iWP = 0; iWP != bTagger.size(); ++iWP) {
    if (debug) cout << "B-Jet: " << bTagger[iWP].first<< " -> " << jet.bDiscriminator(bTagger[iWP].first) <<endl;
   if (jet.bDiscriminator(bTagger[iWP].first) > bTagger[iWP].second ) {
//cout << "Warning, b-tagging called\n";
      taggedJetsHistos[iWP]->analyze(jet, lr, weight);
      if (mcMode) {
	if (!matchB) taggedLightJetsHistos[iWP]->analyze(jet, lr, weight);
	if (matchB) taggedBJetsHistos[iWP]->analyze(jet, lr, weight);
      }
    }
    if (debug) if (!matchB) {
      if (jet.bDiscriminator(bTagger[iWP].first) > bTagger[iWP].second ) cout <<"tagged\n";
      else cout <<"NOT\n";
    }
  }

}
