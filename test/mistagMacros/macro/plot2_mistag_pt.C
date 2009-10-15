#include <fstream.h>
#include <iostream.h>
#include <stdio.h>

plot()
{
// DATA
 TFile *f0 = new TFile("../output/QCD_Pt30-80_7TeV_SD_Jet50U_SSVT.root");
// Primary MC sample
 TFile *f1 = new TFile("../output/QCD_Pt30-80_SD_Jet50U_SSVT.root");
// MC with Track History
 TFile *f2 = new TFile("../output/QCD_Pt80-170_7TeV_SSVT.root");

// *****************************************************************************

 char* xtitle = "p_{T}(jet) (GeV)";
 int nbin = 20;
 float binw = 10., bini = 30.;

 float EtaMin = 0., EtaMax = 2.4;
//  float EtaMin = 0., EtaMax = 0.7;
//  float EtaMin = 0.7, EtaMax = 1.4;
//  float EtaMin = 1.4, EtaMax = 2.4;

//  string TaggerName = "jetProbabilityBJetTags";
//  string TaggerName = "trackCountingHighEffBJetTags";
//  string TaggerName = "trackCountingHighPurBJetTags";
 string TaggerName = "simpleSecondaryVertexBJetTags";

//  float wp = 1.90; // TCHEL
//  float wp = 3.99; // TCHEM
//  float wp = 2.17; // TCHPM
//  float wp = 4.31; // TCHPT
//  float wp = .230; // JPL
//  float wp = .495; // JPM
//  float wp = .700; // JPT
//  float wp = 2.02; // SSVM
 float wp = 3.40; // SSVT

// *****************************************************************************

 int stati=0;
 bool  fit=0;
 bool logy=1;

// *****************************************************************************

 TCanvas *c1 = new TCanvas("c1", "plots",200,10,700,730);
 c1->SetFillColor(10);
 c1->SetFillStyle(4000);
 c1->SetBorderSize(2);

// *****************************************************************************

// TPaveLabel *p01 = new TPaveLabel(0.05,0.93,0.95,0.97,
//                   "mistag in QCD 30-120 CMSSW 1.6.0 : IP2 > 4","br");
// //                   "mistag in QCD 30-120 CMSSW 1.6.0 : IP2 > 4, no pos. IP1 > 4","br");
// // p01->SetFillColor(7);
// p01->SetFillColor(0);
// p01->SetFillStyle(3017);
// p01->SetTextSize(0.8);
// p01->Draw();

pad1 = new TPad("pad1","This is pad1",0.04,0.50,0.96,0.93,21);
pad2 = new TPad("pad2","This is pad2",0.04,0.05,0.96,0.48,21);

pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.05);
   pad1->SetLeftMargin(0.15);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(0);
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.15);
   pad2->SetRightMargin(0.05);
   pad2->SetLeftMargin(0.15);

//$$ gStyle->SetOptDate(1);
gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleColor(29);
gStyle->SetTitleW(0.2);
gStyle->SetTitleH(0.1);
gStyle->SetOptStat(stati);

if (fit) {
gStyle->SetOptFit(111);
gStyle->SetStatW(0.5);
gStyle->SetStatH(0.2);
} else {
gStyle->SetOptFit(0);
gStyle->SetStatW(0.4);
gStyle->SetStatH(0.3);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Systematics

// Flavour fraction
 float Xbb = 0.5;
 float Xcc = 0.5;
 float Xgg = 0.2;

// IP Sign Flip
 float XipM = 0.011 * TMath::Sqrt(1.1);    // 10 pb-1
//  float XipM = 0.011 * TMath::Sqrt(0.2); // 100 pb-1
//  float XipM = 0.011 * TMath::Sqrt(0.11); // 1000 pb-1
 float XipT = 0.014 * TMath::Sqrt(1.1);    // 10 pb-1
//  float XipT = 0.014 * TMath::Sqrt(0.2); // 100 pb-1
//  float XipT = 0.014 * TMath::Sqrt(0.11); // 1000 pb-1
 float Xip = XipM;

// V0 fraction
 float Xv = 0.2;

// Bad tracks fraction
 float Xbad = 0.3;
 float SFbad = 2; // medium

 float BinMin[30], BinMax[30]; 
 float Leff[30], LeffMin[30], LeffMax[30], LeffErr[30];
 float Lsf[30], LsfMin[30], LsfMax[30], LsfErr[30];

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Initialise Histograms

       f1->cd();
 TH1F* g00 = (TH1F*)gROOT->FindObject("hAllFlav_JetPt");   // all
 TH1F* g01 = (TH1F*)gROOT->FindObject("hUDSFlav_JetPt");   // uds-jet
 TH1F* g02 = (TH1F*)gROOT->FindObject("hCFlav_JetPt");     // c-jet
 TH1F* g03 = (TH1F*)gROOT->FindObject("hBFlav_JetPt");     // b-jet
 TH1F* g04 = (TH1F*)gROOT->FindObject("hGluonFlav_JetPt"); // g-jet
 TH1F* g09 = (TH1F*)gROOT->FindObject("hLightFlav_JetPt"); // udsg-jet
//
 TH1F* g10 = (TH1F*)gROOT->FindObject("hAllFlav_NegTag_JetPt");   // all neg
 TH1F* g19 = (TH1F*)gROOT->FindObject("hLightFlav_NegTag_JetPt"); // udsg neg 
//
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Veto_JetPt");   // all with veto
 TH1F* h01 = (TH1F*)gROOT->FindObject("hUDSFlav_Veto_JetPt");   // uds-jet with veto 
 TH1F* h02 = (TH1F*)gROOT->FindObject("hCFlav_Veto_JetPt");     // c-jet with veto 
 TH1F* h03 = (TH1F*)gROOT->FindObject("hBFlav_Veto_JetPt");     // b-jet with veto
 TH1F* h04 = (TH1F*)gROOT->FindObject("hGluonFlav_Veto_JetPt"); // g-jet with veto
 TH1F* h09 = (TH1F*)gROOT->FindObject("hLightFlav_Veto_JetPt"); // udsg-jet with veto 
//
 TH1F* h10 = (TH1F*)gROOT->FindObject("hAllFlav_Veto_NegTag_JetPt");   // all neg with veto
 TH1F* h11 = (TH1F*)gROOT->FindObject("hUDSFlav_Veto_NegTag_JetPt");   // uds neg with veto 
 TH1F* h12 = (TH1F*)gROOT->FindObject("hCFlav_Veto_NegTag_JetPt");     // c neg with veto
 TH1F* h13 = (TH1F*)gROOT->FindObject("hBFlav_Veto_NegTag_JetPt");     // b neg with veto
 TH1F* h14 = (TH1F*)gROOT->FindObject("hGluonFlav_Veto_NegTag_JetPt"); // g neg with veto
//
 TH1F* h20 = (TH1F*)gROOT->FindObject("hAllFlav_PosTag_JetPt");   // all pos
 TH1F* h21 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag_JetPt");   // uds pos 
 TH1F* h22 = (TH1F*)gROOT->FindObject("hCFlav_PosTag_JetPt");     // c pos 
 TH1F* h23 = (TH1F*)gROOT->FindObject("hBFlav_PosTag_JetPt");     // b pos
 TH1F* h24 = (TH1F*)gROOT->FindObject("hGluonFlav_PosTag_JetPt"); // g pos
 TH1F* h29 = (TH1F*)gROOT->FindObject("hLightFlav_PosTag_JetPt"); // udsg pos
//
       g00->Sumw2(); 
       g01->Sumw2(); 
       g02->Sumw2(); 
       g03->Sumw2(); 
       g04->Sumw2(); 
       g09->Sumw2(); 
       g10->Sumw2(); 
       g19->Sumw2(); 
       h00->Sumw2(); 
       h01->Sumw2(); 
       h02->Sumw2(); 
       h03->Sumw2(); 
       h04->Sumw2(); 
       h10->Sumw2(); 
       h11->Sumw2(); 
       h12->Sumw2(); 
       h13->Sumw2(); 
       h14->Sumw2(); 
       h20->Sumw2(); 
       h21->Sumw2(); 
       h22->Sumw2(); 
       h23->Sumw2(); 
       h24->Sumw2(); 
       h29->Sumw2(); 

       f2->cd();
 TH1F* i00 = (TH1F*)gROOT->FindObject("hAllFlav_Veto_JetPt");   // all with veto
 TH1F* i10 = (TH1F*)gROOT->FindObject("hAllFlav_Veto_NegTag_JetPt");   // all neg with veto
 TH1F* i09 = (TH1F*)gROOT->FindObject("hLightFlav_Veto_JetPt"); // udsg-jet with veto 
 TH1F* i29 = (TH1F*)gROOT->FindObject("hLightFlav_PosTag_JetPt"); // udsg pos
//
       i00->Sumw2(); 
       i10->Sumw2(); 
       i09->Sumw2(); 
       i29->Sumw2(); 

 TH1F* h100 = (TH1F*)gROOT->FindObject("hAllFlav_K0s_JetPt");
 TH1F* h109 = (TH1F*)gROOT->FindObject("hLightFlav_K0s_JetPt");
 TH1F* h110 = (TH1F*)gROOT->FindObject("hAllFlav_K0s_Veto_NegTag_JetPt");
 TH1F* h129 = (TH1F*)gROOT->FindObject("hLightFlav_K0s_PosTag_JetPt");
//
 TH1F* h200 = (TH1F*)gROOT->FindObject("hAllFlav_Fak_JetPt");
 TH1F* h209 = (TH1F*)gROOT->FindObject("hLightFlav_Fak_JetPt");
 TH1F* h210 = (TH1F*)gROOT->FindObject("hAllFlav_Fak_Veto_NegTag_JetPt");
 TH1F* h229 = (TH1F*)gROOT->FindObject("hLightFlav_Fak_PosTag_JetPt");
//
 TH1F* h300 = (TH1F*)gROOT->FindObject("hAllFlav_Gam_JetPt");
 TH1F* h309 = (TH1F*)gROOT->FindObject("hLightFlav_Gam_JetPt");
 TH1F* h310 = (TH1F*)gROOT->FindObject("hAllFlav_Gam_Veto_NegTag_JetPt");
 TH1F* h329 = (TH1F*)gROOT->FindObject("hLightFlav_Gam_PosTag_JetPt");
//
       h100->Sumw2(); 
       h109->Sumw2(); 
       h110->Sumw2(); 
       h129->Sumw2(); 
       h200->Sumw2(); 
       h209->Sumw2(); 
       h210->Sumw2(); 
       h229->Sumw2(); 
       h300->Sumw2(); 
       h309->Sumw2(); 
       h310->Sumw2(); 
       h329->Sumw2(); 

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Mistags and Negative Tags

 TH1F* h4= new TH1F("h4","",nbin,bini,bini+nbin*binw);  // light pos / light
       h4->Divide(h29,g09,1,1,"B"); 

 TH1F* h5= new TH1F("h5","",nbin,bini,bini+nbin*binw);  // all neg with veto / all with veto
       h5->Divide(h10,h00,1,1,"B");

 TH1F* Rlight= new TH1F("Rlight","",nbin,bini,bini+nbin*binw);  // (light pos/light) / (all neg/all)
       Rlight->Divide(h4,h5,1,1);

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Flavour fraction

 TH1F* k00= new TH1F("k00","",nbin,bini,bini+nbin*binw);  // g/uds ratio
 TH1F* k01= new TH1F("k01","",nbin,bini,bini+nbin*binw);  // uds-fraction
 TH1F* k02= new TH1F("k02","",nbin,bini,bini+nbin*binw);  // c-fraction
 TH1F* k03= new TH1F("k03","",nbin,bini,bini+nbin*binw);  // b-fraction
 TH1F* k04= new TH1F("k04","flavour fraction",nbin,bini,bini+nbin*binw);  // g-fraction
 TH1F* k11= new TH1F("k11","",nbin,bini,bini+nbin*binw);  // uds-fraction neg
 TH1F* k12= new TH1F("k12","flavour fraction",nbin,bini,bini+nbin*binw);  // c-fraction neg
 TH1F* k13= new TH1F("k13","",nbin,bini,bini+nbin*binw);  // b-fraction neg
 TH1F* k14= new TH1F("k14","",nbin,bini,bini+nbin*binw);  // g-fraction neg
 TH1F* k21= new TH1F("k21","",nbin,bini,bini+nbin*binw);  // uds-fraction pos
 TH1F* k24= new TH1F("k24","",nbin,bini,bini+nbin*binw);  // g-fraction pos

       k00->Divide(g04,g01,1,1);
       k01->Divide(g01,g00,1,1,"B");
       k11->Divide(h11,h10,1,1,"B");
       k21->Divide(h21,h29,1,1,"B");
       k02->Divide(g02,g00,1,1,"B");
       k12->Divide(h12,h10,1,1,"B");
       k03->Divide(g03,g00,1,1,"B");
       k13->Divide(h13,h10,1,1,"B");
       k04->Divide(g04,g00,1,1,"B");
       k14->Divide(h14,h10,1,1,"B");
       k24->Divide(h24,h29,1,1,"B");

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Systematics

// systematics due to flavour fraction

 TH1F* bb= new TH1F("bb","",nbin,bini,bini+nbin*binw);
       for (int i=1; i<nbin+1; i++) {
        float yneg = k13->GetBinContent(i);
	float y = k13->GetBinContent(i) - k03->GetBinContent(i);
	y = Xbb * TMath::Abs( y );
	float ey = Xbb * k13->GetBinError(i);
        bb->SetBinContent(i,y);
        bb->SetBinError(i,ey);
       }
 TH1F* cc= new TH1F("cc","",nbin,bini,bini+nbin*binw);
       for (int i=1; i<nbin+1; i++) {
        float yneg = k12->GetBinContent(i);
	float y = k12->GetBinContent(i) - k02->GetBinContent(i);
	y = Xcc * TMath::Abs( y );
	float ey = Xcc * k12->GetBinError(i);
        cc->SetBinContent(i,y);
        cc->SetBinError(i,ey);
       }

 TH1F* gg= new TH1F("gg","",nbin,bini,bini+nbin*binw);
       for (int i=1; i<nbin+1; i++) {
        float frac = k00->GetBinContent(i);
        float lneg = k11->GetBinContent(i);
        float lpos = k21->GetBinContent(i);
        float gneg = k14->GetBinContent(i);
        float gpos = k24->GetBinContent(i);
	float yl = frac * (lpos - lneg);
	float yg = gpos - gneg;
	y = Xgg * TMath::Abs( yg - yl );
        float elneg = k11->GetBinError(i);
        float elpos = k21->GetBinError(i);
        float egneg = k14->GetBinError(i);
        float egpos = k24->GetBinError(i);
	float eyneg = (egneg - frac*elneg) * (egneg - frac*elneg);
	float eypos = (egpos - frac*elpos) * (egpos - frac*elpos);
	float ey = Xgg * TMath::Sqrt( eyneg + eypos );
        gg->SetBinContent(i,y);
        gg->SetBinError(i,ey);
       }

// systematics due to sign flip

 TH1F* h1= new TH1F("h1","",nbin,bini,bini+nbin*binw); // all pos / all neg
       h1->Divide(h20,g10,1,1);
 TH1F* h2= new TH1F("h2","",nbin,bini,bini+nbin*binw); // light neg / light pos
       h2->Divide(g19,h29,1,1);
 TH1F* sf= new TH1F("sf","",nbin,bini,bini+nbin*binw); // sign flip
       sf->Add(h1,h2,Xip,-Xip);

// systematics due to V0

 TH1F* h71= new TH1F("h71","",nbin,bini,bini+nbin*binw); // K0 neg / all neg
       h71->Divide(h110,i10,1,1,"B");
 TH1F* h72= new TH1F("h72","",nbin,bini,bini+nbin*binw); // light K0 pos / light pos
       h72->Divide(h129,i29,1,1,"B");
 TH1F* h73= new TH1F("h73","",nbin,bini,bini+nbin*binw); // light K0 / light
       h73->Divide(h109,i09,1,1,"B");
 TH1F* h74= new TH1F("h74","",nbin,bini,bini+nbin*binw); // K0 / all
       h74->Divide(h100,i00,1,1,"B");

 TH1F* h81= new TH1F("h81","",nbin,bini,bini+nbin*binw); // Bad neg / all neg
       h81->Divide(h210,i10,SFbad,1,"B");
 TH1F* h82= new TH1F("h82","",nbin,bini,bini+nbin*binw); // light Bad pos / light pos
       h82->Divide(h229,i29,SFbad,1,"B");
 TH1F* h83= new TH1F("h83","",nbin,bini,bini+nbin*binw); // light Bad / light
       h83->Divide(h209,i09,SFbad,1,"B");
 TH1F* h84= new TH1F("h84","",nbin,bini,bini+nbin*binw); // Bad / all
       h84->Divide(h200,i00,SFbad,1,"B");

 TH1F* h91= new TH1F("h91","",nbin,bini,bini+nbin*binw); // gamma neg / all neg
       h91->Divide(h310,i10,1,1,"B");
 TH1F* h92= new TH1F("h92","",nbin,bini,bini+nbin*binw); // light gamma pos / light pos
       h92->Divide(h329,i29,1,1,"B");
 TH1F* h93= new TH1F("h93","",nbin,bini,bini+nbin*binw); // light gamma / light
       h93->Divide(h309,i09,1,1,"B");
 TH1F* h94= new TH1F("h94","",nbin,bini,bini+nbin*binw); // gamma / all
       h94->Divide(h300,i00,1,1,"B");

 TH1F* ka = new TH1F("ka","",nbin,bini,bini+nbin*binw);
 TH1F* bad = new TH1F("bad","",nbin,bini,bini+nbin*binw);
 TH1F* ga = new TH1F("ga","",nbin,bini,bini+nbin*binw);
 TH1F* v0 = new TH1F("v0","",nbin,bini,bini+nbin*binw);

       for (int i=1; i<nbin+1; i++) {
	float y    = h72->GetBinContent(i) - h71->GetBinContent(i);
	             - h73->GetBinContent(i) + h74->GetBinContent(i);
	y = Xv * TMath::Abs( y );
	float ey = h71->GetBinError(i)*h71->GetBinError(i) 
	           + h72->GetBinError(i)*h72->GetBinError(i);
	ey = Xv * TMath::Sqrt( ey );
        ka->SetBinContent(i,y);
        ka->SetBinError(i,ey);
       }
       for (int i=1; i<nbin+1; i++) {
	float y    = h82->GetBinContent(i) - h81->GetBinContent(i);
	             - h83->GetBinContent(i) + h84->GetBinContent(i);
	y = Xbad * TMath::Abs( y );
	float ey = h81->GetBinError(i)*h81->GetBinError(i) 
	           + h82->GetBinError(i)*h82->GetBinError(i);
	ey = Xbad * TMath::Sqrt( ey );
        bad->SetBinContent(i,y);
        bad->SetBinError(i,ey);
       }
       for (int i=1; i<nbin+1; i++) {
	float y    = h92->GetBinContent(i) - h91->GetBinContent(i);
	             - h93->GetBinContent(i) + h94->GetBinContent(i);
	y = Xv * TMath::Abs( y );
	float ey = h91->GetBinError(i)*h91->GetBinError(i) 
	           + h92->GetBinError(i)*h92->GetBinError(i);
	ey = Xv * TMath::Sqrt( ey );
        ga->SetBinContent(i,y);
        ga->SetBinError(i,ey);
       }
       for (int i=1; i<nbin+1; i++) {
	float y = ka->GetBinContent(i)*ka->GetBinContent(i) 
	          + ga->GetBinContent(i)*ga->GetBinContent(i);
	y = TMath::Sqrt( y );
	float ey = ka->GetBinError(i)*ka->GetBinError(i) 
	           + ga->GetBinError(i)*ga->GetBinError(i);
	ey = TMath::Sqrt( ey );
        v0->SetBinContent(i,y);
        v0->SetBinError(i,ey);
       }

 TH1F* tot= new TH1F("tot","",nbin,bini,bini+nbin*binw);
       for (int i=1; i<nbin+1; i++) {
        float ybb = bb->GetBinContent(i);
        float ycc = cc->GetBinContent(i);
        float ygg = gg->GetBinContent(i);
        float ysf = sf->GetBinContent(i);
        float yv0 = v0->GetBinContent(i);
        float ybad = bad->GetBinContent(i);
        float y = ybb + ycc;
	y = TMath::Sqrt( y*y + ygg*ygg + ysf*ysf + yv0*yv0 + ybad*ybad );
        float eybb = bb->GetBinError(i);
        float eycc = cc->GetBinError(i);
        float eygg = gg->GetBinError(i);
        float eysf = sf->GetBinError(i);
        float eyv0 = v0->GetBinError(i);
        float eybad = bad->GetBinError(i);
        float ey = eybb + eycc;
	ey = TMath::Sqrt( ey*ey + eygg*eygg + eysf*eysf + eyv0*eyv0 + eybad*eybad );
        tot->SetBinContent(i,y);
        tot->SetBinError(i,ey);
       }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Mistag in Data = all neg data * Rlight

 TF1 *Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",bini,bini+nbin*binw);
      Fun->SetLineColor(1);
      Fun->SetLineStyle(1);
 TF1 *Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",bini,bini+nbin*binw);
      Gun->SetLineColor(1);
      Gun->SetLineStyle(2);
 TF1 *Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",bini,bini+nbin*binw);
      Hun->SetLineColor(1);
      Hun->SetLineStyle(2);

       f0->cd();
 TH1F* H00 = (TH1F*)gROOT->FindObject("hAllFlav_Veto_JetPt");   // all with veto
 TH1F* H10 = (TH1F*)gROOT->FindObject("hAllFlav_Veto_NegTag_JetPt");   // all neg with veto
       H00->Sumw2(); 
       H10->Sumw2(); 

 TH1F* H5= new TH1F("H5","",nbin,bini,bini+nbin*binw);  // all neg with veto / all with veto
       H5->Divide(H10,H00,1,1,"B");
 TH1F* Mistag= new TH1F("Mistag","",nbin,bini,bini+nbin*binw);  // mistag data
       Mistag->Multiply(H5,Rlight,1,1); 
 TH1F* MisMax= new TH1F("MisMax","",nbin,bini,bini+nbin*binw);  // mistag data up
 TH1F* MisMin= new TH1F("MisMin","",nbin,bini,bini+nbin*binw);  // mistag data down

// MC Truth in "Data"

 TH1F* G09 = (TH1F*)gROOT->FindObject("hLightFlav_JetPt"); // udsg-jet
 TH1F* H29 = (TH1F*)gROOT->FindObject("hLightFlav_PosTag_JetPt"); // udsg pos
       G09->Sumw2(); 
       H29->Sumw2(); 
 TH1F* H4= new TH1F("H4","",nbin,bini,bini+nbin*binw);  // light pos / light
       H4->Divide(H29,G09,1,1,"B"); 

     pad1->cd();
       Mistag->Draw("E"); 
       Mistag->SetLineColor(1);
       Mistag->SetMarkerStyle(20);
       Mistag->SetMarkerColor(1);
       Mistag->SetMarkerSize(1.);
       H4->Draw("Esame");
       H4->SetLineColor(6);
       H4->SetMarkerStyle(21);
       H4->SetMarkerColor(6);
       H4->SetMarkerSize(1.1);
       Mistag->GetXaxis()->SetLabelSize(0.06);
       Mistag->GetYaxis()->SetLabelSize(0.06);
       Mistag->GetXaxis()->SetTitleSize(0.06);
       Mistag->GetXaxis()->SetTitle(xtitle);
       Mistag->GetXaxis()->SetTitleColor(1);
//        Mistag->SetMinimum(0.04); Mistag->SetMaximum(0.4); // Loose
//        Mistag->SetMinimum(0.002); Mistag->SetMaximum(0.05); // Medium
       Mistag->SetMinimum(0.0001); Mistag->SetMaximum(0.02); // Tight
//        Mistag->SetMinimum(0.0001); Mistag->SetMaximum(1.00);
       Mistag->GetXaxis()->SetNdivisions(509);
       Mistag->GetYaxis()->SetNdivisions(509);
       Mistag->Draw("Esame"); 
       Mistag->Fit("Fun","rvee"); 
       H4->Draw("Esame"); 
       Mistag->Draw("Esame"); 
       Fun->Draw("same"); 

       for (int i=1; i<nbin+1; i++) {
	float y = Mistag->GetBinContent(i);
	float ey = Mistag->GetBinError(i);
	float syst = y * tot->GetBinContent(i);
        float stat = Fun->IntegralError(bini+(i-1)*binw,bini+i*binw) / binw;
        float eror = TMath::Sqrt( syst*syst + stat*stat);
	MisMax->SetBinContent(i,y + eror);
	MisMax->SetBinError(i,ey);
	MisMin->SetBinContent(i,y - eror);
	MisMin->SetBinError(i,ey);
	if ( y > 0. ) {
	  MisMax->SetBinError(i,ey * (1. + eror/y ));
	  MisMin->SetBinError(i,ey * (1. - eror/y ));
	}
//$$
    cout << i << " Mistag " << y << " syst " << syst << " stat " << stat << endl;
//$$
	BinMin[i] = bini+(i-1)*binw;
	BinMax[i] = bini+i*binw;
	Leff[i] = Fun->Integral(BinMin[i],BinMax[i]) / binw;;
       }

       MisMax->Fit("Gun","rvee0");
       Gun->Draw("same"); 
       for (int i=1; i<nbin+1; i++) {
	LeffMax[i] = Gun->Integral(BinMin[i],BinMax[i]) / binw;
       }

       MisMin->Fit("Hun","rvee0");
       Hun->Draw("same"); 
       for (int i=1; i<nbin+1; i++) {
	LeffMin[i] = Hun->Integral(BinMin[i],BinMax[i]) / binw;
	LeffErr[i] = (LeffMax[i] - LeffMin[i]) / 2.;
//$$
    cout << BinMin[i] << " " << BinMax[i] << " " 
         << " Light Eff " << int(Leff[i]*1e6)/1e6 << " +_ " << int(LeffErr[i]*1e6)/1e6 << endl;
//$$
       }

  TLegend* leg = new TLegend(0.20,0.73,0.50,0.93);
    leg->SetHeader("Mistag");
    leg->AddEntry(Mistag,"DATA","P");
    leg->AddEntry(H4,"MC truth","P");
    leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Data/MC Scale Factor

 TF1 *FUn = new TF1("FUn","[0]",bini,bini+nbin*binw);
      FUn->SetLineColor(1);
      FUn->SetLineStyle(1);
 TF1 *GUn = new TF1("GUn","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
      GUn->SetLineColor(1);
      GUn->SetLineStyle(2);
 TF1 *HUn = new TF1("HUn","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
      HUn->SetLineColor(1);
      HUn->SetLineStyle(2);

 TH1F* SFlight= new TH1F("SFlight","",nbin,bini,bini+nbin*binw);  // Scale Factor
       SFlight->Divide(H5,h5,1,1); 
 TH1F* SFmax= new TH1F("SFmax","",nbin,bini,bini+nbin*binw);  // SF up
 TH1F* SFmin= new TH1F("SFmin","",nbin,bini,bini+nbin*binw);  // SF down

     pad2->cd();
       SFlight->Draw("E"); 
       SFlight->SetLineColor(1);
       SFlight->SetMarkerStyle(20);
       SFlight->SetMarkerColor(1);
       SFlight->SetMarkerSize(1.);
       SFlight->GetXaxis()->SetLabelSize(0.06);
       SFlight->GetYaxis()->SetLabelSize(0.06);
       SFlight->GetXaxis()->SetTitleSize(0.06);
       SFlight->GetXaxis()->SetTitle(xtitle);
       SFlight->GetXaxis()->SetTitleColor(1);
       SFlight->SetMinimum(0.); 
       SFlight->SetMaximum(2.);
       SFlight->GetXaxis()->SetNdivisions(509);
       SFlight->GetYaxis()->SetNdivisions(509);
       SFlight->Draw("Esame"); 
       SFlight->Fit("FUn","rvee"); 
       FUn->Draw("same"); 

       for (int i=1; i<nbin+1; i++) {
	float y = SFlight->GetBinContent(i);
	float ey = SFlight->GetBinError(i);
	float syst = y * tot->GetBinContent(i);
        float eror = TMath::Sqrt( syst*syst + ey*ey);
        float yfit = FUn->Integral(bini+(i-1)*binw,bini+i*binw) / binw;
	SFmax->SetBinContent(i,yfit + eror);
	SFmax->SetBinError(i,ey);
	SFmin->SetBinContent(i,yfit - eror);
	SFmin->SetBinError(i,ey);
//$$
//$$    cout << i << " SF " << y << " syst " << syst << " stat " << stat << endl;
//$$
	Lsf[i] = FUn->Integral(BinMin[i],BinMax[i]) / binw;
       }

       SFmax->Fit("GUn","rvee0");
       GUn->Draw("same"); 
       for (int i=1; i<nbin+1; i++) {
	LsfMax[i] = GUn->Integral(BinMin[i],BinMax[i]) / binw;
       }

       SFmin->Fit("HUn","rvee0");
       HUn->Draw("same"); 
       for (int i=1; i<nbin+1; i++) {
	LsfMin[i] = HUn->Integral(BinMin[i],BinMax[i]) / binw;
	LsfErr[i] = (LsfMax[i] - LsfMin[i]) / 2.;
//$$
    cout << BinMin[i] << " " << BinMax[i] << " " 
         << " Light SF " << int(Lsf[i]*1e6)/1e6 << " +_ " << int(LsfErr[i]*1e6)/1e6 << endl;
//$$
       }

  TLegend* leg = new TLegend(0.20,0.78,0.50,0.93);
    leg->SetHeader("Scale Factor");
    leg->AddEntry(SFlight,"Neg.Tag Data/MC","P");
    leg->Draw();
    
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Write the ascii file

       string blabla_file = "Mistag.txt";

       ofstream out_txt(blabla_file.c_str(), ios::out);

       out_txt << TaggerName << endl; // < ---the tagger name
       out_txt << wp << endl; // <---- the WP cut on the discriminator
       out_txt << "PerformancePayloadFromTable" << endl;
       out_txt << 4 << endl; // <---- # of results, lighteff and lightefferr
       out_txt << 2 << endl; // <---- # of variables to bin (eta and pt)
       out_txt << "1005 1006 1013 1014" << endl; // <----- results are light eff +_ error, light SF +_ error
       out_txt << "5 2" << endl; // <----- binning is |eta| / pt
       for (int i=1; i<nbin+1; i++) { // <----- eta and Et ranges, light eff +_ error, light SF +_ error
        if (EtaMin == 0 && BinMin[i] < 100 && BinMax[i] < 100) {  
          out_txt << EtaMin    << "       " << EtaMax     << "      " 
                  << BinMin[i] << "      " << BinMax[i]  << "     " 
	          << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
	          << int(Lsf[i]*1e6)/1e6    << "     " << int(LsfErr[i]*1e6)/1e6  << endl;
	}
        else if (EtaMin == 0 && BinMin[i] < 100) {  
          out_txt << EtaMin    << "       " << EtaMax     << "      " 
                  << BinMin[i] << "     " << BinMax[i]  << "     " 
	          << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
	          << int(Lsf[i]*1e6)/1e6    << "     " << int(LsfErr[i]*1e6)/1e6  << endl;
	}
        else if (EtaMin == 0) {  
          out_txt << EtaMin    << "       " << EtaMax     << "     " 
                  << BinMin[i] << "     " << BinMax[i]  << "     " 
	          << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
	          << int(Lsf[i]*1e6)/1e6    << "     " << int(LsfErr[i]*1e6)/1e6  << endl;
	}
        else if (BinMin[i] < 100 && BinMax[i] < 100) {  
          out_txt << EtaMin    << "     " << EtaMax     << "      " 
                  << BinMin[i] << "      " << BinMax[i]  << "     " 
	          << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
	          << int(Lsf[i]*1e6)/1e6    << "     " << int(LsfErr[i]*1e6)/1e6  << endl;
	}
        else if (BinMin[i] < 100) {  
          out_txt << EtaMin    << "     " << EtaMax     << "      " 
                  << BinMin[i] << "     " << BinMax[i]  << "     " 
	          << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
	          << int(Lsf[i]*1e6)/1e6    << "     " << int(LsfErr[i]*1e6)/1e6  << endl;
	}
        else {  
          out_txt << EtaMin    << "     " << EtaMax     << "     " 
                  << BinMin[i] << "     " << BinMax[i]  << "     " 
	          << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
	          << int(Lsf[i]*1e6)/1e6    << "     " << int(LsfErr[i]*1e6)/1e6  << endl;
	}
       }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  c1->Update();
}

