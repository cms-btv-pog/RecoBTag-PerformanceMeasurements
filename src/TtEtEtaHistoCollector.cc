#include "RecoBTag/PerformanceMeasurements/interface/TtEtEtaHistoCollector.h"
// #include "RecoBTag/Analysis/interface/JetTagPlotter.h"
// #include "FWCore/Utilities/interface/CodedException.h"
#include "TStyle.h"

using namespace reco;
using namespace std;

TtEtEtaHistoCollector::TtEtEtaHistoCollector(const TString & baseName, 
	const std::vector<double> & etRanges, 
	const std::vector<double> & etaRanges, const double lowerBound,
	const double upperBound, const int bins, const bool update) :
	etRanges_(etRanges), etaRanges_(etaRanges), lowerBound_(lowerBound),
	upperBound_(upperBound), bins_(bins)
{
  baseName_ = baseName;
  diffPlotOK = false;

  if (update){
    gFile->cd();
    gFile->cd(baseName_);
  }

  //
  // Book all histograms.
  //

  // iterate over ranges:
  const int iEtaStart = -1                   ;  // this will be the inactive one
  const int iEtaEnd   = etaRanges.size() - 1 ;
  const int iEtStart  = -1                   ;  // this will be the inactive one
  const int iEtEnd    = etRanges.size() - 1  ;
//   setTDRStyle();

    // the objects for the differential plots vs. eta,et for

    vector<TtDifferentialPlot*> * differentialPlotsConstantEta = new vector<TtDifferentialPlot*> () ;
    for ( int iEta = iEtaStart ; iEta < iEtaEnd ; iEta++ ) {
      TtDifferentialPlot * etaConstTtDifferentialPlot = new TtDifferentialPlot
	(this, TtDifferentialPlot::constETA, baseName);
      differentialPlotsConstantEta->push_back ( etaConstTtDifferentialPlot );
    }

    vector<TtDifferentialPlot*> * differentialPlotsConstantEt  = new vector<TtDifferentialPlot*> () ;
    for ( int iEt = iEtStart ; iEt < iEtEnd ; iEt++ ) {
      // differentialPlots for this et bin
      TtDifferentialPlot * etConstTtDifferentialPlot = new TtDifferentialPlot
	(this, TtDifferentialPlot::constET, baseName);
      differentialPlotsConstantEt->push_back ( etConstTtDifferentialPlot );
    }

    // eta loop
  for ( int iEta = iEtaStart ; iEta < iEtaEnd ; iEta++ ) {
    // et loop
    for ( int iEt = iEtStart ; iEt < iEtEnd ; iEt++ ) {

      // DEFINE BTagBin:
      bool    etaActive_ , etActive_;
      double  etaMin_, etaMax_, etMin_, etMax_ ;

      if ( iEta != -1 ) {
	etaActive_ = true ;
	etaMin_    = etaRanges[iEta]   ;
	etaMax_    = etaRanges[iEta+1] ;
      }
      else {
	etaActive_ = false ;
	etaMin_    = etaRanges[0]   ;
	etaMax_    = etaRanges[etaRanges.size() - 1]   ;
      }

      if ( iEt != -1 ) {
	etActive_ = true ;
	etMin_    = etRanges[iEt]   ;
	etMax_    = etRanges[iEt+1] ;
      }
      else {
	etActive_ = false ;
	etMin_    = etRanges[0]	;
	etMax_    = etRanges[etRanges.size() - 1]	;
      }
      EtaEtBin * etaEtBin = new EtaEtBin(etaActive_ , etaMin_ , etaMax_ ,
			etActive_  , etMin_  , etMax_ );

      // Instantiate the genertic b tag plotter

      TString name = baseName + etaEtBin->getDescriptionString();
      TH1F* histo;
      if (!update) {
	histo = new TH1F ( name, name, bins , lowerBound , upperBound );
	histo->Sumw2() ; 
      } else {
        histo   = (TH1F*) gDirectory->Get(name) ; 
      }
      binHistoMap[etaEtBin] = histo;

	// Add to the corresponding differential plotters
	(*differentialPlotsConstantEta)[iEta+1]->addBin ( etaEtBin ) ;
	(*differentialPlotsConstantEt )[iEt+1] ->addBin ( etaEtBin ) ;
    }
  }

      // the objects for the differential plots vs. eta, et: collect all from constant eta and constant et
      differentialPlots.reserve(differentialPlotsConstantEta->size()+differentialPlotsConstantEt->size()) ;
      differentialPlots.insert(differentialPlots.end(), differentialPlotsConstantEta->begin(), differentialPlotsConstantEta->end());
      differentialPlots.insert(differentialPlots.end(), differentialPlotsConstantEt->begin(), differentialPlotsConstantEt->end());

//       edm::LogInfo("Info")
// 	cout << "====>>>> ## sizeof differentialPlots = " << differentialPlots.size();

      // the intermediate ones are no longer needed
      delete differentialPlotsConstantEta ;
      delete differentialPlotsConstantEt  ;

}



// TtEtEtaHistoCollector::TtEtEtaHistoCollector(const TtEtEtaHistoCollector & other)
// {
//   for (BinHistoMapConstIter bin = other.histoBegin();
// 	bin!=other.histoEnd(); ++bin) {
//       binHistoMap[new EtaEtBin(*bin->first)] = (TH1F*)bin->second->Clone("");
//   }
//   for (vector<double>::const_iterator bin = other.etRanges_.begin();
// 	bin!=other.etRanges_.end(); ++bin) {
//       etRanges_.push_back(*bin);
//   }
//   for (vector<double>::const_iterator bin = other.etaRanges_.begin();
// 	bin!=other.etaRanges_.end(); ++bin) {
//       etaRanges_.push_back(*bin);
//   }
// 
//   diffPlotVect differentialPlots;
// 
//   lowerBound_ = other.lowerBound_;
//   upperBound_ = other.upperBound_;
//   bins_ = other.bins_;
// 
//   baseName_ = other.baseName();
//   
// //   copy diff
//   
// }
// 
// TtEtEtaHistoCollector & TtEtEtaHistoCollector::operator=
// 	(TtEtEtaHistoCollector & other)
// {
//   etRanges_.clear();
//   etaRanges_.clear();
//   for (BinHistoMap::iterator iBin = binHistoMap.begin();
// 	iBin != binHistoMap.end(); ++iBin) {
//     delete iBin->second;
//   }
//   binHistoMap.clear();
//   for (BinHistoMapIter bin = other.histoBegin();
// 	bin!=other.histoEnd(); ++bin) {
//       binHistoMap[new EtaEtBin(*bin->first)] = (TH1F*)bin->second->Clone("");
//   }
//   for (vector<double>::iterator bin = other.etRanges_.begin();
// 	bin!=other.etRanges_.end(); ++bin) {
//       etRanges_.push_back(*bin);
//   }
//   for (vector<double>::iterator bin = other.etaRanges_.begin();
// 	bin!=other.etaRanges_.end(); ++bin) {
//       etaRanges_.push_back(*bin);
//   }
//   lowerBound_ = other.lowerBound_;
//   upperBound_ = other.upperBound_;
//   bins_ = other.bins_;
// 
// //   copy diff
// 
// 
//   baseName_ = other.baseName();
//   return *this;
// }


TtEtEtaHistoCollector::~TtEtEtaHistoCollector()
{
  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
    delete iBin->second;
  }
  for (vector<TtDifferentialPlot *>::iterator iPlotter = differentialPlots.begin();
    	      iPlotter != differentialPlots.end(); ++ iPlotter) {
    delete *iPlotter;
  }
}


void TtEtEtaHistoCollector::Sumw2()
{
  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
    iBin->second->Sumw2();
  }
}

void TtEtEtaHistoCollector::analyze(const TopJet & jet, const double lhr,
	const double weight)
{
  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
    if ( iBin->first->inBin(jet) ) 
      iBin->second->Fill(lhr, weight);
  }
}


void TtEtEtaHistoCollector::write()
{
  gFile->cd();
  gFile->mkdir(baseName_);
  gFile->cd(baseName_);
    cout << "write :"<<baseName_<<endl;

//   setTDRStyle();
  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
      iBin->second->Write();
//       if (producePs)  (*binPlotters[iPlotter]).psPlot(psBaseName);
//       if (produceEps) (*binPlotters[iPlotter]).epsPlot(epsBaseName);
    }

  for (vector<TtDifferentialPlot *>::iterator iPlotter = differentialPlots.begin();
    	      iPlotter != differentialPlots.end(); ++ iPlotter) {
//      (**iPlotter).process();
//       if (producePs)  (**iPlotter).psPlot(psBaseName);
//       if (produceEps) (**iPlotter).epsPlot(epsBaseName);
      (**iPlotter).write();
  }
  gFile->cd();
}


void TtEtEtaHistoCollector::buildTtDifferentialPlots(int bin)
{
  for (vector<TtDifferentialPlot *>::iterator iPlotter = differentialPlots.begin();
    	      iPlotter != differentialPlots.end(); ++ iPlotter) {
     (**iPlotter).process(bin);
  }
  diffPlotOK = true;
}

void TtEtEtaHistoCollector::setHisto(const EtaEtBin& bin, TH1F * histo)
{
  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
    if (*(iBin->first) == bin) {
      delete iBin->second;
      iBin->second = histo;
    }
  }
}

TH1F* TtEtEtaHistoCollector::getHisto(const EtaEtBin& bin)
{
  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
    if (*(iBin->first) == bin) return iBin->second;
  }
  return 0;
}

TtEtEtaHistoCollector * TtEtEtaHistoCollector::buildIntegratedCollector()
{
  TtEtEtaHistoCollector * newColl = new TtEtEtaHistoCollector("Int_"+baseName_ , etRanges_, etaRanges_,
	lowerBound_, upperBound_, bins_, false);
    for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
      integratedHisto(iBin->second, newColl->getHisto(*(iBin->first)) );
    }
  return newColl;
}

TtEtEtaHistoCollector * TtEtEtaHistoCollector::buildRatioCollector(const TString& name,
			TtEtEtaHistoCollector & otherColl)
{
  TtEtEtaHistoCollector * newColl = new TtEtEtaHistoCollector(name, etRanges_, etaRanges_,
	lowerBound_, upperBound_, bins_, false);
    for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
      newColl->getHisto(*(iBin->first))->Divide(iBin->second, 
	otherColl.getHisto(*(iBin->first)), 1., 1., "B");
    }
  return newColl;
}

TtEtEtaHistoCollector * TtEtEtaHistoCollector::buildDiffCollector(const TString& name,
			TtEtEtaHistoCollector & otherColl)
{
  TtEtEtaHistoCollector * newColl = new TtEtEtaHistoCollector(name, etRanges_, etaRanges_,
	lowerBound_, upperBound_, bins_, false);
    for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
      newColl->getHisto(*(iBin->first))->Add(iBin->second, 
	otherColl.getHisto(*(iBin->first)), -1.);
    }
  return newColl;
}


// TH1F* TtEtEtaHistoCollector::integratedHisto(const EtaEtBin& bin)
// {
//   for (BinHistoMap::iterator iBin = binHistoMap.begin();
// 	iBin != binHistoMap.end(); ++iBin) {
//     if (*(iBin->first) == bin) {
//       TString name = baseName_ +"_Int_"+ iBin->first->getDescriptionString();
//       return integratedHisto(iBin->second, name);
//     }
//   }
//   return 0;
// }

void TtEtEtaHistoCollector::integratedHisto(const TH1F * const histo, TH1F * newHisto)
{
  int bins = histo->GetNbinsX();
  double integral;
//   TH1F* intHisto = new TH1F ( newName, newName, bins, histo->GetBinLowEdge(0),
//    histo->GetBinLowEdge(bins+1) );
  for(int bin=0; bin<=bins+1; ++bin){
    integral = histo->Integral(bin,bins+1);
    newHisto->SetBinContent(bin, integral);
    newHisto->SetBinError(bin, sqrt(integral) );
  //fLRtotSoverSplusB->Eval(LRcutVal);
//    float Sint = hLRtotS->Integral(cut, hLRtotS->GetNbinsX()+1);
//    float Bint = hLRtotB->Integral(cut, hLRtotB->GetNbinsX()+1);
//    cout << cut << " " <<Eff[cut] << " " <<hLRtotS->Integral(0,cut)
//    <<" " <<hLRtotSoverSplusB->GetBinContent(cut)
//    <<" " << Sint <<" " << Bint<<" " << Sint/(Sint+Bint)
//    <<endl;
  }
}

void TtEtEtaHistoCollector::eps(TH1F * histo, const TString & ext)
{
  TStyle *tdrStyle = gROOT->GetStyle("tdrStyle");
  tdrStyle->SetOptFit(0);
  tdrStyle->SetOptStat(0);
  tdrStyle->SetOptTitle(0);
//   tdrStyle->SetPadTopMargin(0.01);
//   tdrStyle->SetPadBottomMargin(0.01);
//   tdrStyle->SetPadLeftMargin(0.01);
//   tdrStyle->SetPadRightMargin(0.01);

  TCanvas c2("c2","",600,600);
//   c2.SetTopMargin(0.01);
//   c2.SetBottomMargin(0.01);
//   c2.SetLeftMargin(0.01);
//   c2.SetRightMargin(0.01);
  histo->Draw("P");
  cout << histo->GetName()+TString(".")+ext<<endl;
  c2.Print(histo->GetName()+TString(".")+ext);
}

void TtEtEtaHistoCollector::eps(const TString & ext)
{
  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
    eps(iBin->second, ext);
  }
  if (diffPlotOK) 
    for (vector<TtDifferentialPlot *>::iterator iPlotter = differentialPlots.begin();
    	      iPlotter != differentialPlots.end(); ++ iPlotter) {
      (**iPlotter).eps(ext);
  }
}

void TtEtEtaHistoCollector::setRange(const float min, const float max)
{

  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
   iBin->second->GetYaxis()->SetRangeUser(min, max);
  }
  if (diffPlotOK) 
    for (vector<TtDifferentialPlot *>::iterator iPlotter = differentialPlots.begin();
    	      iPlotter != differentialPlots.end(); ++ iPlotter) {
      (**iPlotter).setRange(min, max);
  }
}


void TtEtEtaHistoCollector::setLabels(const TString& xLabel, const TString& yLabel)
{

  for (BinHistoMap::iterator iBin = binHistoMap.begin();
	iBin != binHistoMap.end(); ++iBin) {
    iBin->second->GetXaxis()->SetTitle(xLabel);
    iBin->second->GetYaxis()->SetTitle(yLabel);
  }
  if (diffPlotOK) 
    for (vector<TtDifferentialPlot *>::iterator iPlotter = differentialPlots.begin();
    	      iPlotter != differentialPlots.end(); ++ iPlotter) {
      (**iPlotter).setLabels(yLabel);
  }
}
