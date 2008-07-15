#include <iostream>
#include <cassert>
#include <TROOT.h>
#include <TSystem.h>
#include <Cintex/Cintex.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TKey.h>
#include <vector>
using namespace std;

TString  fileName    = "events.root";


int main( int argc, const char* argv[]) { 
  double W=0, Z=0, tt=0, QCD=0;

  // TString fileName = dir+dataFile;
  bool verbose = false;
  for (int i=1;i<argc;++i) {
    if (strncmp(argv[i],"-v",2)==0) verbose=true;
    else fileName = argv[i];
  }

  cout <<"File : "<<fileName<<endl;
  TFile *theFile = new TFile (fileName) ;
  TH1F* histo   = (TH1F*) gDirectory->Get("ProcessesAfter") ; 
  if (verbose) {
    cout << "Events per channel:\n";
    for (int i=1;i<histo->GetNbinsX(); i++) {
      cout <<i<<" "<< histo->GetXaxis()->GetBinLabel(i) <<"\t"<<histo->GetBinContent(i)<<endl;
    }
  }
  
  for (int i=24;i<=28; i++) tt+=histo->GetBinContent(i);
  cout << "Top Events (24-28): \t"<<tt<<endl;
  for (int i=2;i<=12; i++) W+=histo->GetBinContent(i);
  cout << "W   Events ( 2-12): \t"<<W<<endl;
  for (int i=13;i<=23; i++) Z+=histo->GetBinContent(i);
  cout << "Z   Events (13-23): \t"<<Z<<endl;
  for (int i=30;i<=49; i++) QCD+=histo->GetBinContent(i);
  QCD+=histo->GetBinContent(63);
  cout << "QCD Events (30-49 & 63): \t"<<QCD<<endl;
  cout <<"Other\t"<<histo->Integral(1,71)-tt-W-Z-QCD<<endl;
  theFile->Close();
}
