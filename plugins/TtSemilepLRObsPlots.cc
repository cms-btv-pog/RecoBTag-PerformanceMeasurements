#include "RecoBTag/PerformanceMeasurements/plugins/TtSemilepLRObsPlots.h"
#include "CSA07EffAnalyser/CSA07EffAnalyser/interface/CSA07ProcessId.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"

using namespace std;
using namespace reco;
using namespace math;

//
// constructors and destructor
//
TtSemilepLRObsPlots::TtSemilepLRObsPlots(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  debug = iConfig.getParameter<bool> ("debug");
  rootFileName      = iConfig.getParameter<string> ("rootFileName"); 
  csa = iConfig.getParameter< bool > ("CSA");
  if (!csa) {
    weight = iConfig.getParameter< double > ("weight");
  }
 
  nrJetCombObs   = iConfig.getParameter<int> ("nrJetCombObs");
  obsNrs	   = iConfig.getParameter< vector<int> > ("JetCombObs");
  nrJetCombHistBins	   = iConfig.getParameter< int > ("nrJetCombHistBins");
  obsMin	   = iConfig.getParameter< vector<double> > ("JetCombObsMin");
  obsMax	   = iConfig.getParameter< vector<double> > ("JetCombObsMax");

  evtsols    = iConfig.getParameter<edm::InputTag>("evtSolution");
  jetSource_ = iConfig.getParameter<edm::InputTag>("jetSource");
  
  bTagCutLabel = iConfig.getParameter< string >("bTagCutLabel");
  bCut = iConfig.getParameter< double >("bCut");

  for(int j = 0; j < nrJetCombObs; j++){
    obsFits.push_back("pol4");
  }
  
  theFile = new TFile("a.root", "RECREATE");
  theFile->cd();

  // define all histograms & fit functions
  myLRhelper = new LRHelpFunctions(obsNrs, nrJetCombHistBins, obsMin, obsMax, obsFits);  

  allSolution = 0;
  goodSolution2 = 0;
  goodSolution1 = 0;
  goodSolution = 0;
  B=0;nonB=0;tau=0;

  myLRJetCombObservables = new TtSemiLRJetCombObservables();
  myLRJetCombObservables->jetSource(jetSource_);
}


TtSemilepLRObsPlots::~TtSemilepLRObsPlots()
{
   delete myLRhelper;
   theFile->Close();
}

void TtSemilepLRObsPlots::endJob()
{
  cout << "TtSemiLepLRObsPlots: Events per fb-1 " << goodSolution2 << " w/o btagcut w/o chi2cut"<< endl;
  cout << "TtSemiLepLRObsPlots: Events per fb-1 " << goodSolution1 << " btagcut w/o chi2cut"<< endl;
  cout << "TtSemiLepLRObsPlots: Events per fb-1 " << goodSolution << " btagcut & chi2cut"<< endl;
  //    pair<double, double> match = myLRhelper->getBnumbers();
  // cout << "TtSemiLepLRObsPlots match " << B*weight<<" "<<nonB*weight<< " "<< tau<<" "<<match.second<<endl;
  
  // store histograms and fits in root-file
  myLRhelper -> storeToROOTfile(rootFileName);
  if (debug) cout << "************* Finished writing histograms to file" << endl;
  
}

void TtSemilepLRObsPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get the event solution
  edm::Handle< std::vector<TtSemiEvtSolution> > eSols; 
  iEvent.getByLabel(evtsols, eSols);
  
  std::vector<TtSemiEvtSolution> sols;
  sols = *eSols;
  if (debug) cout << "TtSemiLepLRObsPlots: Found "<< sols.size()  << " Semilepton solutions.\n";
 
  if(sols.size()== 12) {

   if (csa) {  
   if (debug) cout << "CSA07 sample" << endl;
     //reading out the event weights needed for processing the CSA07 soups
      edm::Handle< double> weightHandle;
      iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
      weight = * weightHandle;

      //it would be nice to have a detailed counting of the differen jet multiplicity bins and different channels

      //int procID = csa07::csa07ProcessId(iEvent);
      //if (debug) cout << "processID: " << procID
      //<< " - name: " << csa07::csa07ProcessName(procID) 
      //<< " - weight: "<< weight << endl;
    }

    double maxProbChi2=-999;
    double maxbTagCut=-9999;

    for (int w = 0; w < 12; w++){  
      //extra selection cuts
      if (sols[w].getProbChi2()>maxProbChi2){maxProbChi2=sols[w].getProbChi2();}   
      if (sols[w].getCalHadb().getBDiscriminator(bTagCutLabel) > maxbTagCut) {maxbTagCut=sols[w].getCalHadb().getBDiscriminator(bTagCutLabel);} 
    }
    if (debug) cout << maxProbChi2 << " " << maxbTagCut << endl;
    
    goodSolution2 += weight;
    if(maxbTagCut>bCut)  {
      goodSolution1 += weight;
      if (maxProbChi2>0){
        goodSolution += weight;
        if(debug) cout << "this event is selected" <<endl;
      }
    }
    
   // if(maxProbChi2>0&&maxbTagCut>bCut)  {
   //   goodSolution += weight;
   //   if(debug) cout << "this event is selected" <<endl;
   // }
    
    int matchMode =0;
    vector < vector< TtSemiLRJetCombObservables::IntBoolPair > > obsMatch;
    vector <int> matchSum;
    int bestSol=-1, bestMatch = -1;
    int corrWeight=0;
    vector <int> goodMatch(nrJetCombObs); 
    vector <int> badMatch(nrJetCombObs); 
    vector <int> totMatch(nrJetCombObs);
    
    for (int i = 0; i < 12; ++i) {
      obsMatch.push_back(myLRJetCombObservables->operator()(sols[i],iEvent, true));
    }	  
    
    if (debug) cout << "--- b-Tag cut --- " << endl;
    if (debug) cout << bTagCutLabel << " cut at : " << bCut << endl;
    
    for(int j = 0; j < nrJetCombObs; ++j) { 
      goodMatch[j]=0;
      badMatch[j]=0;
      totMatch[j]=0;
    }
    
    if (debug) cout << "Start calculating weights for observables depending on the matching" << endl;
    // Fill the observables 
    if(matchMode==0){
      for(int j = 0; j < nrJetCombObs; ++j) { 
	for (int w = 0; w < 12; w++){      //this 2 selection cuts can only be applied solution by solution 
	  //if (debug) cout << "Chi2 prob :  " << sols[w].getProbChi2() << endl;
	  if (sols[w].getProbChi2()>0) { // don't consider non converged solutions
	    //if (debug) cout << "b-Tag value : " << sols[w].getCalHadb().getBDiscriminator(bTagCutLabel) << endl;
	    if (sols[w].getCalHadb().getBDiscriminator(bTagCutLabel) > bCut) { // b-tag cut to suppres W+jets
	      if( myLRhelper->obsFitIncluded(obsNrs[j]) ) { 
		if (debug) cout << "Matching : " << obsMatch[w][obsNrs[j]-1].second<<endl;
		goodMatch[j]+=obsMatch[w][obsNrs[j]-1].second;
		badMatch[j]+=!obsMatch[w][obsNrs[j]-1].second;
		totMatch[j]++;
	      }  
	    }
	  }
	}
      }
    }
    
    
    
    
    /*  if (debug) {
	for(int j = 0; j < nrJetCombObs; ++j) { 
	if( myLRhelper->obsFitIncluded(obsNrs[j]) ) { 
	cout << "obs " << j << ": goodmatch " << goodMatch[j] << "  badmatch " << badMatch[j] << "  totmatch " << totMatch[j] << endl;
	}  
	}
	}*/
    
    if (matchMode==1){
      for(int j = 0; j < nrJetCombObs; ++j) { 
	for (int w = 0; w < 12; w++){      //this 2 selection cuts can only be applied solution by solution 
	  if (debug) cout << "Chi2 prob :  " << sols[w].getProbChi2() << endl;
	  if (sols[w].getProbChi2()>0) { // don't consider non converged solutions
	    if (debug) cout << "b-Tag value : " << sols[w].getCalHadb().getBDiscriminator(bTagCutLabel) << endl;
	    if (sols[w].getCalHadb().getBDiscriminator(bTagCutLabel) > bCut) { // b-tag cut to suppres W+jets
	      if( myLRhelper->obsFitIncluded(obsNrs[j]) ) { 
		if (fabs(sols[w].getCalLepb().getPartonFlavour())==5) {goodMatch[j]++;} else {badMatch[j]++;}
		totMatch[j]++;
	      }  
	    }
	  }
	}
      }
    }
    
    
    for (int i = 0; i < 12; ++i) { 
      
      //this 2 selection cuts can only be applied solution by solution
      if (sols[i].getProbChi2()>0) { // don't consider non converged solutions
	if (sols[i].getCalHadb().getBDiscriminator(bTagCutLabel) > bCut) { // b-tag cut to suppres W+jets
	  
	  // Fill the observables 
	  for(int j = 0; j < nrJetCombObs; ++j) {
	    if( myLRhelper->obsFitIncluded(obsNrs[j]) ) {
	      
	      if (matchMode == 0){
		//FIXME: dirty hack to make obs 62 (ET4/ET5) not normalized
		bool sandbPresent = false;
		if (goodMatch[j]!=0&&badMatch[j]!=0&&j!=61){ sandbPresent = true;}
		if(j==61){sandbPresent = true;}
		
		if(sandbPresent){
		  if (obsMatch[i][obsNrs[j]-1].second) {
		    
		    myLRhelper->fillToSignalHists(obsNrs[j],sols[i].getLRJetCombObsVal(obsNrs[j]), weight*totMatch[j]/goodMatch[j]);
		    
		    for(int k = j; k < nrJetCombObs; ++k) {
		      if (obsMatch[i][obsNrs[k]-1].second) {
			
			for (int w = 0; w < 12; w++){
			  corrWeight+=(obsMatch[w][obsNrs[j]-1].second * obsMatch[w][obsNrs[k]-1].second);
			} 
			
			myLRhelper->fillToSignalCorrelation(obsNrs[j], sols[i].getLRJetCombObsVal(obsNrs[j]), obsNrs[k], sols[i].getLRJetCombObsVal(obsNrs[k]), weight/corrWeight); //FIXME: probably corrWeight wrongly calculated
			
			// cout << "Fill :" << obsMatch[0][obsNrs[j]-1].second*obsMatch[0][obsNrs[k]-1].second<<" "<<
			// 		  	obsMatch[0][obsNrs[j]-1].second<<" "<<obsMatch[0][obsNrs[k]-1].second<<" "<<
			// 			obsNrs[j]<<" "<< sols[i].getLRSignalEvtObsVal(obsNrs[j])<<" "<<
			// obsNrs[k]<<" "<< sols[i].getLRSignalEvtObsVal(obsNrs[k])<<endl;
		      }
		    }
		  } else {
		  myLRhelper->fillToBackgroundHists(obsNrs[j],sols[i].getLRJetCombObsVal(obsNrs[j]), weight*totMatch[j]/badMatch[j]);
		  }
		}
	      }//end of MatchMode
	      
	      //***new matchMode
	      if (matchMode == 1){
		//if (goodMatch[j]!=0&&badMatch[j]!=0){
		if (sols[i].getCalLepb().getPartonFlavour()==5) {
		  
		  myLRhelper->fillToSignalHists(obsNrs[j],sols[i].getLRJetCombObsVal(obsNrs[j]), weight*totMatch[j]/goodMatch[j]);
		  
		  for(int k = j; k < nrJetCombObs; ++k) {
		    if (sols[i].getCalLepb().getPartonFlavour()==5) {
		      
		      for (int w = 0; w < 12; w++){
			//corrWeight+=(obsMatch[w][obsNrs[j]-1].second * obsMatch[w][obsNrs[k]-1].second);
			corrWeight+=1;//FIXME
		      } 
		      
		      myLRhelper->fillToSignalCorrelation(obsNrs[j], sols[i].getLRJetCombObsVal(obsNrs[j]), obsNrs[k], sols[i].getLRJetCombObsVal(obsNrs[k]), weight/corrWeight); //FIXME: probably corrWeight wrongly calculated
		      
		      // cout << "Fill :" << obsMatch[0][obsNrs[j]-1].second*obsMatch[0][obsNrs[k]-1].second<<" "<<
		      // 		  	obsMatch[0][obsNrs[j]-1].second<<" "<<obsMatch[0][obsNrs[k]-1].second<<" "<<
		      // 			obsNrs[j]<<" "<< sols[i].getLRSignalEvtObsVal(obsNrs[j])<<" "<<
		      // obsNrs[k]<<" "<< sols[i].getLRSignalEvtObsVal(obsNrs[k])<<endl;
		    }
		  }
		} else myLRhelper->fillToBackgroundHists(obsNrs[j],sols[i].getLRJetCombObsVal(obsNrs[j]), weight*totMatch[j]/badMatch[j]);
		//	}
	      }//end of MatchMode
	      
	      
	    }
	  }
	  
	  
	  int matchSum = 0;
	  for(int j = 0; j < nrJetCombObs; j++){
	    if (obsMatch[i][j].second) ++matchSum;
	    if (debug) cout <<obsMatch[i][j].second;
	  }
	  if (debug) cout <<endl;
	  if (matchSum>bestMatch) {
	    bestMatch = matchSum;
	    bestSol = i;
	  }
	}
      }
    }
    //FIXME: fix this cout for 12 solution case and for the extra cuts I've made
    /*if (debug) cout << "\nBest Solution:" <<sols[0].getLRBestJetComb() << " - "<<sols[1].getLRBestJetComb()<< " = "<<bestMatch<<" = "<<bestSol<< endl;
      
    for(int i = 0; i < 12; i++){
    B+=obsMatch[bestSol][i].second;
    nonB+=(!obsMatch[bestSol][i].second);
    }
    if (debug)  cout<< B <<" " <<nonB << " " << obsMatch[bestSol][0].second << " " <<obsMatch[bestSol][1].second<<endl;*/
    //FIXME
    //if (sols[bestSol].getJetB().getPartonFlavour()==5) ++tau;
    //if (sols[bestSol].getJetBbar().getPartonFlavour()==5) ++tau;
    
  }
  if (debug) cout <<" ============================ End TtSemiLepLRObsPlots ============================" <<endl;
}


bool TtSemilepLRObsPlots::allMatch(const TtSemiEvtSolution & sol) 
{/*
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
 */ 
  return true;
 
}

