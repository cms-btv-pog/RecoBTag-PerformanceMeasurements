#include <iostream>
#include <TROOT.h>

plot()
{
// DATA
 TFile *f0 = new TFile("../output/Data_SSVHPT.root");
// Primary MC sample
 TFile *f1 = new TFile("../output/MC_QCD-Pt153080_SSVHPT.root");
// MC with Track History
 TFile *f2 = new TFile("../output/MC_QCD-Pt153080_SSVHPT.root");

// *****************************************************************************

 char* htitle = "SSVHPT"; 
 char* xtitle = "p_{T}(jet) (GeV)";
 int nbin = 10;
 float binw = 10., bini = 10.;

// *****************************************************************************

 Int_t stati=0;
 Bool_t  fit=0;
 Bool_t logy=0;

// *****************************************************************************

 TCanvas *c1 = new TCanvas("c1", "plots",200,10,700,750);
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

pad1 = new TPad("pad1","This is pad1",0.04,0.49,0.48,0.91,21);
pad2 = new TPad("pad2","This is pad2",0.52,0.49,0.96,0.91,21);
pad3 = new TPad("pad3","This is pad5",0.04,0.04,0.48,0.46,21);
pad4 = new TPad("pad4","This is pad6",0.52,0.04,0.96,0.46,21);

pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(0);
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.02);
   pad1->SetLeftMargin(0.18);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(0);
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.15);
   pad2->SetRightMargin(0.02);
   pad2->SetLeftMargin(0.18);

pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(logy);
   pad3->SetTopMargin(0.05);
   pad3->SetBottomMargin(0.15);
   pad3->SetRightMargin(0.02);
   pad3->SetLeftMargin(0.18);

pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(0);
   pad4->SetTopMargin(0.05);
   pad4->SetBottomMargin(0.15);
   pad4->SetRightMargin(0.02);
   pad4->SetLeftMargin(0.18);

//$$ gStyle->SetOptDate(1);
gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(42);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.05);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
// gStyle->SetTitleX(0); // Set the position of the title box
// gStyle->SetTitleY(0.985); // Set the position of the title box
// gStyle->SetTitleStyle(Style_t style = 1001);
// gStyle->SetTitleBorderSize(2);
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
 float Xbb = 0.2;
 float Xcc = 0.5;
 float Xgg = 0.2;

// // IP Sign Flip
//  float XipM = 0.011 * TMath::Sqrt(1.1);    // 10 pb-1
// //  float XipM = 0.011 * TMath::Sqrt(0.2); // 100 pb-1
// //  float XipM = 0.011 * TMath::Sqrt(0.11); // 1000 pb-1
//  float XipT = 0.014 * TMath::Sqrt(1.1);    // 10 pb-1
// //  float XipT = 0.014 * TMath::Sqrt(0.2); // 100 pb-1
// //  float XipT = 0.014 * TMath::Sqrt(0.11); // 1000 pb-1
//  float Xip = XipM;

// V0 fraction
//  float Xk0 = 0.25;
//  float SFk0 = 1.95;
//  float Xla = 0.4;
//  float SFla = 2.2;
 float Xv = 0.25;
 float SFv = 1.95;

// Bad tracks uncertainty
 float Xbad = 0.5;

 float BinMin[30], BinMax[30]; 
 float Leff[30], LeffMin[30], LeffMax[30], LeffErr[30];
 float Lsf[30], LsfMin[30], LsfMax[30], LsfErr[30];

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Initialise Histograms

       f1->cd();
//
 TH1F* h00 = (TH1F*)gROOT->FindObject("hData_JetPt");      // all
 TH1F* h01 = (TH1F*)gROOT->FindObject("hUDSFlav_JetPt");   // uds-jet
 TH1F* h02 = (TH1F*)gROOT->FindObject("hCFlav_JetPt");     // c-jet
 TH1F* h03 = (TH1F*)gROOT->FindObject("hBFlav_JetPt");     // b-jet
 TH1F* h04 = (TH1F*)gROOT->FindObject("hGluonFlav_JetPt"); // g-jet
 TH1F* h09 = (TH1F*)gROOT->FindObject("hLightFlav_JetPt"); // udsg-jet
//
 TH1F* h10 = (TH1F*)gROOT->FindObject("hData_NegTag_JetPt");      // all neg
 TH1F* h11 = (TH1F*)gROOT->FindObject("hUDSFlav_NegTag_JetPt");   // uds neg
 TH1F* h12 = (TH1F*)gROOT->FindObject("hCFlav_NegTag_JetPt");    // c neg
 TH1F* h13 = (TH1F*)gROOT->FindObject("hBFlav_NegTag_JetPt");     // b neg
 TH1F* h14 = (TH1F*)gROOT->FindObject("hGluonFlav_NegTag_JetPt"); // g neg
 TH1F* h19 = (TH1F*)gROOT->FindObject("hLightFlav_NegTag_JetPt"); // udsg neg 
//
 TH1F* h20 = (TH1F*)gROOT->FindObject("hData_PosTag_JetPt");      // all pos
 TH1F* h21 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag_JetPt");   // uds pos 
 TH1F* h22 = (TH1F*)gROOT->FindObject("hCFlav_PosTag_JetPt");     // c pos 
 TH1F* h23 = (TH1F*)gROOT->FindObject("hBFlav_PosTag_JetPt");     // b pos
 TH1F* h24 = (TH1F*)gROOT->FindObject("hGluonFlav_PosTag_JetPt"); // g pos
 TH1F* h29 = (TH1F*)gROOT->FindObject("hLightFlav_PosTag_JetPt"); // udsg pos
//
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
       h19->Sumw2(); 
       h20->Sumw2(); 
       h21->Sumw2(); 
       h22->Sumw2(); 
       h23->Sumw2(); 
       h24->Sumw2(); 
       h29->Sumw2(); 

       f2->cd();
 TH1F* i00 = (TH1F*)gROOT->FindObject("hData_JetPt");             // all
 TH1F* i10 = (TH1F*)gROOT->FindObject("hData_NegTag_JetPt");      // all neg
 TH1F* i09 = (TH1F*)gROOT->FindObject("hLightFlav_JetPt");        // udsg-jet
 TH1F* i29 = (TH1F*)gROOT->FindObject("hLightFlav_PosTag_JetPt"); // udsg pos
       i00->Sumw2(); 
       i10->Sumw2(); 
       i09->Sumw2(); 
       i29->Sumw2(); 

 TH1F* h100 = (TH1F*)gROOT->FindObject("hAllFlav_K0s_JetPt");
 TH1F* h109 = (TH1F*)gROOT->FindObject("hLightFlav_K0s_JetPt");
 TH1F* h110 = (TH1F*)gROOT->FindObject("hAllFlav_K0s_NegTag_JetPt");
 TH1F* h129 = (TH1F*)gROOT->FindObject("hLightFlav_K0s_PosTag_JetPt");
//
 TH1F* h200 = (TH1F*)gROOT->FindObject("hAllFlav_Fak_JetPt");
 TH1F* h209 = (TH1F*)gROOT->FindObject("hLightFlav_Fak_JetPt");
 TH1F* h210 = (TH1F*)gROOT->FindObject("hAllFlav_Fak_NegTag_JetPt");
 TH1F* h229 = (TH1F*)gROOT->FindObject("hLightFlav_Fak_PosTag_JetPt");
//
 TH1F* h300 = (TH1F*)gROOT->FindObject("hAllFlav_Gam_JetPt");
 TH1F* h309 = (TH1F*)gROOT->FindObject("hLightFlav_Gam_JetPt");
 TH1F* h310 = (TH1F*)gROOT->FindObject("hAllFlav_Gam_NegTag_JetPt");
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
       h4->Divide(h29,h09,1,1,"B"); 
 TH1F* h5= new TH1F("h5","",nbin,bini,bini+nbin*binw);  // all neg / all
       h5->Divide(h10,h00,1,1,"B");
 TH1F* Rlight= new TH1F("Rlight",htitle,nbin,bini,bini+nbin*binw);  // (light pos/light) / (all neg/all)
       Rlight->Divide(h4,h5,1,1);

 TH1F* H4 = new TH1F("H4", "",nbin,bini,bini+nbin*binw);
 TH1F* H5 = new TH1F("H5", "",nbin,bini,bini+nbin*binw);
 TH1F* H00= new TH1F("H00","",nbin,bini,bini+nbin*binw);
 TH1F* H10= new TH1F("H10","",nbin,bini,bini+nbin*binw);
 TH1F* H29= new TH1F("H29","",nbin,bini,bini+nbin*binw);
 TH1F* H09= new TH1F("H09","",nbin,bini,bini+nbin*binw);
       H00->Sumw2(); 
       H10->Sumw2(); 
       H29->Sumw2(); 
       H09->Sumw2(); 
       H00->Add(h00,h100,1,SFv-1); 
       H10->Add(h10,h110,1,SFv-1); 
       H29->Add(h29,h129,1,SFv-1); 
       H09->Add(h09,h109,1,SFv-1); 
       H4->Divide(H29,H09,1,1,"B"); 
       H5->Divide(H10,H00,1,1,"B");
 TH1F* RlightWeight= new TH1F("RlightWeight",htitle,nbin,bini,bini+nbin*binw);  // (light pos/light) / (all neg/all)
       RlightWeight->Divide(H4,H5,1,1);

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

       k00->Divide(h04,h01,1,1);
       k01->Divide(h01,h00,1,1,"B");
       k11->Divide(h11,h10,1,1,"B");
       k21->Divide(h21,h29,1,1,"B");
       k02->Divide(h02,h00,1,1,"B");
       k12->Divide(h12,h10,1,1,"B");
       k03->Divide(h03,h00,1,1,"B");
       k13->Divide(h13,h10,1,1,"B");
       k04->Divide(h04,h00,1,1,"B");
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

// // systematics due to sign flip
// 
//  TH1F* h1= new TH1F("h1","",nbin,bini,bini+nbin*binw); // all pos / all neg
//        h1->Divide(h20,h10,1,1);
//  TH1F* h2= new TH1F("h2","",nbin,bini,bini+nbin*binw); // light neg / light pos
//        h2->Divide(h19,h29,1,1);
//  TH1F* sf= new TH1F("sf","",nbin,bini,bini+nbin*binw); // sign flip
//        sf->Add(h1,h2,Xip,-Xip);

// systematics due to K0s+Lambda
 TH1F* h71= new TH1F("h71","",nbin,bini,bini+nbin*binw); // K0 neg / all neg
       h71->Divide(h110,i10,1,1,"B");
 TH1F* h72= new TH1F("h72","",nbin,bini,bini+nbin*binw); // light K0 pos / light pos
       h72->Divide(h129,i29,1,1,"B");
 TH1F* h73= new TH1F("h73","",nbin,bini,bini+nbin*binw); // light K0 / light
       h73->Divide(h109,i09,1,1,"B");
 TH1F* h74= new TH1F("h74","",nbin,bini,bini+nbin*binw); // K0 / all
       h74->Divide(h100,i00,1,1,"B");

// systematics due to fakes
float n8 = hLightFlav_PosTag_Fak->GetEntries();
float n10= hLightFlav_PosTag_Bad->GetEntries() - hLightFlav_NegTag_Bad->GetEntries();
    n10 += hLightFlav_PosTag_Oth->GetEntries() - hLightFlav_NegTag_Oth->GetEntries();
 float SFbad = 1 + n10 / n8;
 TH1F* h81= new TH1F("h81","",nbin,bini,bini+nbin*binw); // Bad neg / all neg
       h81->Divide(h210,i10,SFbad,1,"B");
 TH1F* h82= new TH1F("h82","",nbin,bini,bini+nbin*binw); // light Bad pos / light pos
       h82->Divide(h229,i29,SFbad,1,"B");
 TH1F* h83= new TH1F("h83","",nbin,bini,bini+nbin*binw); // light Bad / light
       h83->Divide(h209,i09,SFbad,1,"B");
 TH1F* h84= new TH1F("h84","",nbin,bini,bini+nbin*binw); // Bad / all
       h84->Divide(h200,i00,SFbad,1,"B");

// systematics due to 2ndary interactions
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
//         float ysf = sf->GetBinContent(i);
        float yv0 = v0->GetBinContent(i);
        float ybad = bad->GetBinContent(i);
        float y = ybb + ycc;
// 	y = TMath::Sqrt( y*y + ygg*ygg + ysf*ysf + yv0*yv0 + ybad*ybad );
	y = TMath::Sqrt( y*y + ygg*ygg + yv0*yv0 + ybad*ybad );
        float eybb = bb->GetBinError(i);
        float eycc = cc->GetBinError(i);
        float eygg = gg->GetBinError(i);
//         float eysf = sf->GetBinError(i);
        float eyv0 = v0->GetBinError(i);
        float eybad = bad->GetBinError(i);
        float ey = eybb + eycc;
// 	ey = TMath::Sqrt( ey*ey + eygg*eygg + eysf*eysf + eyv0*eyv0 + eybad*eybad );
	ey = TMath::Sqrt( ey*ey + eygg*eygg + eyv0*eyv0 + eybad*eybad );
        tot->SetBinContent(i,y);
        tot->SetBinError(i,ey);
       }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Mistag in Data = all neg data * Rlight

//  TF1 *Fun = new TF1("Fun","[0]+[1]*x",bini,bini+nbin*binw);
 TF1 *Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      Fun->SetLineColor(1);
      Fun->SetLineStyle(1);
      Fun->SetLineWidth(0.5);
// TF1 *Gun = new TF1("Gun","[0]+[1]*x",bini,bini+nbin*binw);
  TF1 *Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      Gun->SetLineColor(1);
      Gun->SetLineStyle(2);
//  TF1 *Hun = new TF1("Hun","[0]+[1]*x",bini,bini+nbin*binw);
 TF1 *Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      Hun->SetLineColor(1);
      Hun->SetLineStyle(2);

       f0->cd();
 TH1F* H00 = (TH1F*)gROOT->FindObject("hData_JetPt");	     // all
 TH1F* H10 = (TH1F*)gROOT->FindObject("hData_NegTag_JetPt"); // all neg
       H00->Sumw2(); 
       H10->Sumw2(); 

 TH1F* Negtag= new TH1F("Negtag","",nbin,bini,bini+nbin*binw); // all neg / all
       Negtag->Divide(H10,H00,1,1,"B");
 TH1F* Mistag= new TH1F("Mistag","",nbin,bini,bini+nbin*binw); // mistag data
       Mistag->Multiply(Negtag,RlightWeight,1,1); 

 TH1F* Mis   = new TH1F("Mis"   ,"",nbin/2,bini,bini+nbin*binw);
 TH1F* MisMax= new TH1F("MisMax","",nbin/2,bini,bini+nbin*binw); // mistag data up
 TH1F* MisMin= new TH1F("MisMin","",nbin/2,bini,bini+nbin*binw); // mistag data down
       for (int i=2; i<nbin+1; i = i+2) {
        float y  = Mistag->GetBinContent(i-1) + Mistag->GetBinContent(i);
        float ey = Mistag->GetBinError(i-1)*Mistag->GetBinError(i-1)
	         + Mistag->GetBinError(i)*Mistag->GetBinError(i);
	ey = TMath::Sqrt(ey);
	int j = i/2;
        Mis->SetBinContent(j,y);
        Mis->SetBinError(j,ey);
       }

     pad1->cd();
       Mis->Draw("E"); 
       Mis->SetLineColor(2);
       Mis->SetMarkerStyle(20);
       Mis->SetMarkerColor(2);
       Mis->SetMarkerSize(1.);
       Mis->SetMinimum(0.0); Mis->SetMaximum(0.010); // tight
       Mis->SetTitleFont(42, "XYZ");
       Mis->GetXaxis()->SetTitle(xtitle);
       Mis->GetYaxis()->SetTitle("mistag rate");
       Mis->GetYaxis()->SetTitleOffset(1.6);
       Mis->SetTickLength(0.03, "XYZ");
       Mis->SetLabelOffset(0.007,"XYZ");
       Mis->SetLabelSize(0.06, "XYZ");
       Mis->SetLabelFont(42, "XYZ");
       Mis->SetTitleSize(0.06, "XYZ");
       Mis->SetNdivisions(509,"XYZ");
       Mis->Draw("Esame"); 
       Mis->Fit("Fun","rvee"); 
       Mis->Draw("Esame"); 
       Fun->Draw("same"); 

       for (int i=1; i<nbin/2+1; i++) {
	float y = Mis->GetBinContent(i);
	float ey = Mis->GetBinError(i);
	float syst = y * tot->GetBinContent(i);
        float stat = Fun->IntegralError(bini+(i-1)*binw*2,bini+i*binw*2) /binw/2;
        float eror = TMath::Sqrt( syst*syst + stat*stat);
	MisMax->SetBinContent(i,y + eror);
	MisMax->SetBinError(i,ey);
	MisMin->SetBinContent(i,y - eror);
	MisMin->SetBinError(i,ey);
	if ( y > 0. ) {
	  MisMax->SetBinError(i,ey * (1. + eror/y ));
	  MisMin->SetBinError(i,ey * (1. - eror/y ));
	}
	BinMin[i] = bini+(i-1)*binw*2;
	BinMax[i] = bini+i*binw*2;
	Leff[i] = Fun->Integral(BinMin[i],BinMax[i]) /binw/2;;
       }

       MisMax->Fit("Gun","rvee0");
       Gun->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LeffMax[i] = Gun->Integral(BinMin[i],BinMax[i]) / binw/2;
       }

       MisMin->Fit("Hun","rvee0");
       Hun->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LeffMin[i] = Hun->Integral(BinMin[i],BinMax[i]) / binw/2;
	LeffErr[i] = (LeffMax[i] - LeffMin[i]) / 2.;
       }

  //include text
  TPaveText *abc = new TPaveText(0.22,0.69,0.60,0.89,"brNDC");   
  abc->SetBorderSize(0);	
  abc->SetFillStyle(0);  
  abc->SetTextAlign(13); 
  abc->SetTextFont(42);	
  abc->SetTextSize(0.07);   
  TText *text = abc->AddText("SSVHPT tagger");   
  TText *text = abc->AddText("|#eta(jet)| < 2.4");   
  abc->Draw(""); 

  //include the official CMS label
  TPaveText *pt = new TPaveText(0.49,0.94,0.99,0.99,"brNDC");   
  pt->SetBorderSize(0);   
  pt->SetFillStyle(0);  
  pt->SetTextAlign(13);   
  pt->SetTextFont(42);   
  pt->SetTextSize(0.05);   
  TText *text = pt->AddText("CMS Preliminary 2010");   
  pt->Draw(""); 

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Data/MC Scale Factor

 TF1 *FUn = new TF1("FUn","([0]+[1]*x)",bini,bini+nbin*binw);
//  TF1 *FUn = new TF1("FUn","([0]+[1]*x)/(1+[2]*x)",bini,bini+nbin*binw);
//  TF1 *FUn = new TF1("FUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      FUn->SetLineColor(1);
      FUn->SetLineStyle(1);
      FUn->SetLineWidth(0.5);
 TF1 *GUn = new TF1("GUn","([0]+[1]*x)",bini,bini+nbin*binw);
//  TF1 *GUn = new TF1("GUn","([0]+[1]*x)/(1+[2]*x)",bini,bini+nbin*binw);
//  TF1 *GUn = new TF1("GUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      GUn->SetLineColor(1);
      GUn->SetLineStyle(2);
 TF1 *HUn = new TF1("HUn","([0]+[1]*x)",bini,bini+nbin*binw);
//  TF1 *HUn = new TF1("HUn","([0]+[1]*x)/(1+[2]*x)",bini,bini+nbin*binw);
//  TF1 *HUn = new TF1("HUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      HUn->SetLineColor(1);
      HUn->SetLineStyle(2);

 TH1F* SFlight= new TH1F("SFlight","",nbin,bini,bini+nbin*binw);  // Scale Factor
       SFlight->Divide(Mistag,h4,1,1); 
 TH1F* SFneg= new TH1F("SFneg","",nbin,bini,bini+nbin*binw);  // neg tag Data/MC
       SFneg->Divide(Negtag,h5,1,1); 
       for (int i=1; i<nbin+1; i++) {
	float ey = SFneg->GetBinError(i);
	SFlight->SetBinError(i,ey);
       }

 TH1F* SF    = new TH1F("SF"  ,"",nbin/2,bini,bini+nbin*binw);
 TH1F* SFmax= new TH1F("SFmax","",nbin/2,bini,bini+nbin*binw);  // SF up
 TH1F* SFmin= new TH1F("SFmin","",nbin/2,bini,bini+nbin*binw);  // SF down
       for (int i=2; i<nbin+1; i = i+2) {
        float y  = SFlight->GetBinContent(i-1) + SFlight->GetBinContent(i);
        float ey = SFlight->GetBinError(i-1)*SFlight->GetBinError(i-1)
	         + SFlight->GetBinError(i)*SFlight->GetBinError(i);
	y = y / 2.;
	ey = TMath::Sqrt(ey) / 2.;
	int j = i/2;
        SF->SetBinContent(j,y);
        SF->SetBinError(j,ey);
       }

     pad2->cd();
       SF->Draw("E"); 
       SF->SetLineColor(2);
       SF->SetMarkerStyle(20);
       SF->SetMarkerColor(2);
       SF->SetMarkerSize(1.);
       SF->SetMinimum(0.0); 
       SF->SetMaximum(3.0);
       SF->SetTitleFont(42, "XYZ");
       SF->GetXaxis()->SetTitle(xtitle);
       SF->GetYaxis()->SetTitle("Data/MC");
       SF->GetYaxis()->SetTitleOffset(1.2);
       SF->SetTickLength(0.03, "XYZ");
       SF->SetLabelOffset(0.007,"XYZ");
       SF->SetLabelSize(0.06, "XYZ");
       SF->SetLabelFont(42, "XYZ");
       SF->SetTitleSize(0.06, "XYZ");
       SF->SetNdivisions(509,"XYZ");
       SF->Draw("Esame"); 
       SF->Fit("FUn","rvee"); 
       FUn->Draw("same"); 
 
       for (int i=1; i<nbin/2+1; i++) {
	float y = SF->GetBinContent(i);
	float ey = SF->GetBinError(i);
	float syst = y * tot->GetBinContent(i);
        float eror = TMath::Sqrt( syst*syst + ey*ey);
        float yfit = FUn->Integral(bini+(i-1)*binw*2,bini+i*binw*2) / binw/2;
	SFmax->SetBinContent(i,yfit + eror);
	SFmax->SetBinError(i,ey);
	SFmin->SetBinContent(i,yfit - eror);
	SFmin->SetBinError(i,ey);
	Lsf[i] = FUn->Integral(BinMin[i],BinMax[i]) / binw/2;
       }

       SFmax->Fit("GUn","rvee0");
       GUn->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LsfMax[i] = GUn->Integral(BinMin[i],BinMax[i]) / binw/2;
       }

       SFmin->Fit("HUn","rvee0");
       HUn->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LsfMin[i] = HUn->Integral(BinMin[i],BinMax[i]) / binw/2;
	LsfErr[i] = (LsfMax[i] - LsfMin[i]) / 2.;
       }

  //include text
  TPaveText *abc = new TPaveText(0.22,0.69,0.60,0.89,"brNDC");   
  abc->SetBorderSize(0);	
  abc->SetFillStyle(0);  
  abc->SetTextAlign(13); 
  abc->SetTextFont(42);	
  abc->SetTextSize(0.07);   
  TText *text = abc->AddText("SSVHPT tagger");   
  TText *text = abc->AddText("|#eta(jet)| < 2.4");   
  abc->Draw(""); 

  //include the official CMS label
  TPaveText *pt = new TPaveText(0.49,0.94,0.99,0.99,"brNDC");   
  pt->SetBorderSize(0);   
  pt->SetFillStyle(0);  
  pt->SetTextAlign(13);   
  pt->SetTextFont(42);   
  pt->SetTextSize(0.05);   
  TText *text = pt->AddText("CMS Preliminary 2010");   
  pt->Draw(""); 

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// DATA
 TFile *f0 = new TFile("../output/Data_ptGT30_SSVHPT.root");
// Primary MC sample
 TFile *f1 = new TFile("../output/MC_QCD-Pt15_ptGT30_SSVHPT.root");
// MC with Track History
 TFile *f2 = new TFile("../output/MC_QCD-Pt15_ptGT30_SSVHPT.root");

// *****************************************************************************

 char* htitle = "SSVHPT"; 
 char* xtitle = "|#eta(jet)|";
 int nbin = 12;
 float binw = 0.2, bini = 0.;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Initialise Histograms

       f1->cd();
//
 TH1F* h00 = (TH1F*)gROOT->FindObject("hData_JetEta");      // all
 TH1F* h01 = (TH1F*)gROOT->FindObject("hUDSFlav_JetEta");   // uds-jet
 TH1F* h02 = (TH1F*)gROOT->FindObject("hCFlav_JetEta");     // c-jet
 TH1F* h03 = (TH1F*)gROOT->FindObject("hBFlav_JetEta");     // b-jet
 TH1F* h04 = (TH1F*)gROOT->FindObject("hGluonFlav_JetEta"); // g-jet
 TH1F* h09 = (TH1F*)gROOT->FindObject("hLightFlav_JetEta"); // udsg-jet
//
 TH1F* h10 = (TH1F*)gROOT->FindObject("hData_NegTag_JetEta");      // all neg
 TH1F* h11 = (TH1F*)gROOT->FindObject("hUDSFlav_NegTag_JetEta");   // uds neg
 TH1F* h12 = (TH1F*)gROOT->FindObject("hCFlav_NegTag_JetEta");    // c neg
 TH1F* h13 = (TH1F*)gROOT->FindObject("hBFlav_NegTag_JetEta");     // b neg
 TH1F* h14 = (TH1F*)gROOT->FindObject("hGluonFlav_NegTag_JetEta"); // g neg
 TH1F* h19 = (TH1F*)gROOT->FindObject("hLightFlav_NegTag_JetEta"); // udsg neg 
//
 TH1F* h20 = (TH1F*)gROOT->FindObject("hData_PosTag_JetEta");      // all pos
 TH1F* h21 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag_JetEta");   // uds pos 
 TH1F* h22 = (TH1F*)gROOT->FindObject("hCFlav_PosTag_JetEta");     // c pos 
 TH1F* h23 = (TH1F*)gROOT->FindObject("hBFlav_PosTag_JetEta");     // b pos
 TH1F* h24 = (TH1F*)gROOT->FindObject("hGluonFlav_PosTag_JetEta"); // g pos
 TH1F* h29 = (TH1F*)gROOT->FindObject("hLightFlav_PosTag_JetEta"); // udsg pos
//
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
       h19->Sumw2(); 
       h20->Sumw2(); 
       h21->Sumw2(); 
       h22->Sumw2(); 
       h23->Sumw2(); 
       h24->Sumw2(); 
       h29->Sumw2(); 

       f2->cd();
 TH1F* i00 = (TH1F*)gROOT->FindObject("hData_JetEta");             // all
 TH1F* i10 = (TH1F*)gROOT->FindObject("hData_NegTag_JetEta");      // all neg
 TH1F* i09 = (TH1F*)gROOT->FindObject("hLightFlav_JetEta");        // udsg-jet
 TH1F* i29 = (TH1F*)gROOT->FindObject("hLightFlav_PosTag_JetEta"); // udsg pos
       i00->Sumw2(); 
       i10->Sumw2(); 
       i09->Sumw2(); 
       i29->Sumw2(); 

 TH1F* h100 = (TH1F*)gROOT->FindObject("hAllFlav_K0s_JetEta");
 TH1F* h109 = (TH1F*)gROOT->FindObject("hLightFlav_K0s_JetEta");
 TH1F* h110 = (TH1F*)gROOT->FindObject("hAllFlav_K0s_NegTag_JetEta");
 TH1F* h129 = (TH1F*)gROOT->FindObject("hLightFlav_K0s_PosTag_JetEta");
//
 TH1F* h200 = (TH1F*)gROOT->FindObject("hAllFlav_Fak_JetEta");
 TH1F* h209 = (TH1F*)gROOT->FindObject("hLightFlav_Fak_JetEta");
 TH1F* h210 = (TH1F*)gROOT->FindObject("hAllFlav_Fak_NegTag_JetEta");
 TH1F* h229 = (TH1F*)gROOT->FindObject("hLightFlav_Fak_PosTag_JetEta");
//
 TH1F* h300 = (TH1F*)gROOT->FindObject("hAllFlav_Gam_JetEta");
 TH1F* h309 = (TH1F*)gROOT->FindObject("hLightFlav_Gam_JetEta");
 TH1F* h310 = (TH1F*)gROOT->FindObject("hAllFlav_Gam_NegTag_JetEta");
 TH1F* h329 = (TH1F*)gROOT->FindObject("hLightFlav_Gam_PosTag_JetEta");
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
       h4->Divide(h29,h09,1,1,"B"); 
 TH1F* h5= new TH1F("h5","",nbin,bini,bini+nbin*binw);  // all neg / all
       h5->Divide(h10,h00,1,1,"B");
 TH1F* Rlight= new TH1F("Rlight",htitle,nbin,bini,bini+nbin*binw);  // (light pos/light) / (all neg/all)
       Rlight->Divide(h4,h5,1,1);

 TH1F* H4 = new TH1F("H4", "",nbin,bini,bini+nbin*binw);
 TH1F* H5 = new TH1F("H5", "",nbin,bini,bini+nbin*binw);
 TH1F* H00= new TH1F("H00","",nbin,bini,bini+nbin*binw);
 TH1F* H10= new TH1F("H10","",nbin,bini,bini+nbin*binw);
 TH1F* H29= new TH1F("H29","",nbin,bini,bini+nbin*binw);
 TH1F* H09= new TH1F("H09","",nbin,bini,bini+nbin*binw);
       H00->Sumw2(); 
       H10->Sumw2(); 
       H29->Sumw2(); 
       H09->Sumw2(); 
       H00->Add(h00,h100,1,SFv-1); 
       H10->Add(h10,h110,1,SFv-1); 
       H29->Add(h29,h129,1,SFv-1); 
       H09->Add(h09,h109,1,SFv-1); 
       H4->Divide(H29,H09,1,1,"B"); 
       H5->Divide(H10,H00,1,1,"B");
 TH1F* RlightWeight= new TH1F("RlightWeight",htitle,nbin,bini,bini+nbin*binw);  // (light pos/light) / (all neg/all)
       RlightWeight->Divide(H4,H5,1,1);

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

       k00->Divide(h04,h01,1,1);
       k01->Divide(h01,h00,1,1,"B");
       k11->Divide(h11,h10,1,1,"B");
       k21->Divide(h21,h29,1,1,"B");
       k02->Divide(h02,h00,1,1,"B");
       k12->Divide(h12,h10,1,1,"B");
       k03->Divide(h03,h00,1,1,"B");
       k13->Divide(h13,h10,1,1,"B");
       k04->Divide(h04,h00,1,1,"B");
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

// // systematics due to sign flip
// 
//  TH1F* h1= new TH1F("h1","",nbin,bini,bini+nbin*binw); // all pos / all neg
//        h1->Divide(h20,h10,1,1);
//  TH1F* h2= new TH1F("h2","",nbin,bini,bini+nbin*binw); // light neg / light pos
//        h2->Divide(h19,h29,1,1);
//  TH1F* sf= new TH1F("sf","",nbin,bini,bini+nbin*binw); // sign flip
//        sf->Add(h1,h2,Xip,-Xip);

// systematics due to K0s+Lambda
 TH1F* h71= new TH1F("h71","",nbin,bini,bini+nbin*binw); // K0 neg / all neg
       h71->Divide(h110,i10,1,1,"B");
 TH1F* h72= new TH1F("h72","",nbin,bini,bini+nbin*binw); // light K0 pos / light pos
       h72->Divide(h129,i29,1,1,"B");
 TH1F* h73= new TH1F("h73","",nbin,bini,bini+nbin*binw); // light K0 / light
       h73->Divide(h109,i09,1,1,"B");
 TH1F* h74= new TH1F("h74","",nbin,bini,bini+nbin*binw); // K0 / all
       h74->Divide(h100,i00,1,1,"B");

// systematics due to fakes
float n8 = hLightFlav_PosTag_Fak->GetEntries();
float n10= hLightFlav_PosTag_Bad->GetEntries() - hLightFlav_NegTag_Bad->GetEntries();
    n10 += hLightFlav_PosTag_Oth->GetEntries() - hLightFlav_NegTag_Oth->GetEntries();
 float SFbad = 1 + n10 / n8;
 TH1F* h81= new TH1F("h81","",nbin,bini,bini+nbin*binw); // Bad neg / all neg
       h81->Divide(h210,i10,SFbad,1,"B");
 TH1F* h82= new TH1F("h82","",nbin,bini,bini+nbin*binw); // light Bad pos / light pos
       h82->Divide(h229,i29,SFbad,1,"B");
 TH1F* h83= new TH1F("h83","",nbin,bini,bini+nbin*binw); // light Bad / light
       h83->Divide(h209,i09,SFbad,1,"B");
 TH1F* h84= new TH1F("h84","",nbin,bini,bini+nbin*binw); // Bad / all
       h84->Divide(h200,i00,SFbad,1,"B");

// systematics due to 2ndary interactions
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
//         float ysf = sf->GetBinContent(i);
        float yv0 = v0->GetBinContent(i);
        float ybad = bad->GetBinContent(i);
        float y = ybb + ycc;
// 	y = TMath::Sqrt( y*y + ygg*ygg + ysf*ysf + yv0*yv0 + ybad*ybad );
	y = TMath::Sqrt( y*y + ygg*ygg + yv0*yv0 + ybad*ybad );
        float eybb = bb->GetBinError(i);
        float eycc = cc->GetBinError(i);
        float eygg = gg->GetBinError(i);
//         float eysf = sf->GetBinError(i);
        float eyv0 = v0->GetBinError(i);
        float eybad = bad->GetBinError(i);
        float ey = eybb + eycc;
// 	ey = TMath::Sqrt( ey*ey + eygg*eygg + eysf*eysf + eyv0*eyv0 + eybad*eybad );
	ey = TMath::Sqrt( ey*ey + eygg*eygg + eyv0*eyv0 + eybad*eybad );
        tot->SetBinContent(i,y);
        tot->SetBinError(i,ey);
       }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Mistag in Data = all neg data * Rlight

//  TF1 *Fun = new TF1("Fun","[0]+[1]*x",bini,bini+nbin*binw);
 TF1 *Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      Fun->SetLineColor(1);
      Fun->SetLineStyle(1);
      Fun->SetLineWidth(0.5);
//  TF1 *Gun = new TF1("Gun","[0]+[1]*x",bini,bini+nbin*binw);
 TF1 *Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      Gun->SetLineColor(1);
      Gun->SetLineStyle(2);
//  TF1 *Hun = new TF1("Hun","[0]+[1]*x",bini,bini+nbin*binw);
 TF1 *Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x",bini,bini+nbin*binw);
//  TF1 *Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      Hun->SetLineColor(1);
      Hun->SetLineStyle(2);

       f0->cd();
 TH1F* H00 = (TH1F*)gROOT->FindObject("hData_JetEta");	     // all
 TH1F* H10 = (TH1F*)gROOT->FindObject("hData_NegTag_JetEta"); // all neg
       H00->Sumw2(); 
       H10->Sumw2(); 

 TH1F* Negtag= new TH1F("Negtag","",nbin,bini,bini+nbin*binw); // all neg / all
       Negtag->Divide(H10,H00,1,1,"B");
 TH1F* Mistag= new TH1F("Mistag","",nbin,bini,bini+nbin*binw); // mistag data
       Mistag->Multiply(Negtag,RlightWeight,1,1); 

 TH1F* Mis   = new TH1F("Mis"   ,"",nbin/2,bini,bini+nbin*binw);
 TH1F* MisMax= new TH1F("MisMax","",nbin/2,bini,bini+nbin*binw); // mistag data up
 TH1F* MisMin= new TH1F("MisMin","",nbin/2,bini,bini+nbin*binw); // mistag data down
       for (int i=2; i<nbin+1; i = i+2) {
        float y  = Mistag->GetBinContent(i-1) + Mistag->GetBinContent(i);
        float ey = Mistag->GetBinError(i-1)*Mistag->GetBinError(i-1)
	         + Mistag->GetBinError(i)*Mistag->GetBinError(i);
	ey = TMath::Sqrt(ey);
	int j = i/2;
        Mis->SetBinContent(j,y);
        Mis->SetBinError(j,ey);
       }

     pad3->cd();
       Mis->Draw("E"); 
       Mis->SetLineColor(2);
       Mis->SetMarkerStyle(20);
       Mis->SetMarkerColor(2);
       Mis->SetMarkerSize(1.);
       Mis->SetMinimum(0.0); Mis->SetMaximum(0.010); // tight
       Mis->SetTitleFont(42, "XYZ");
       Mis->GetXaxis()->SetTitle(xtitle);
       Mis->GetYaxis()->SetTitle("mistag rate");
       Mis->GetYaxis()->SetTitleOffset(1.6);
       Mis->SetTickLength(0.03, "XYZ");
       Mis->SetLabelOffset(0.007,"XYZ");
       Mis->SetLabelSize(0.06, "XYZ");
       Mis->SetLabelFont(42, "XYZ");
       Mis->SetTitleSize(0.06, "XYZ");
       Mis->SetNdivisions(509,"XYZ");
       Mis->Draw("Esame"); 
       Mis->Fit("Fun","rvee"); 
       Mis->Draw("Esame"); 
       Fun->Draw("same"); 

       for (int i=1; i<nbin/2+1; i++) {
	float y = Mis->GetBinContent(i);
	float ey = Mis->GetBinError(i);
	float syst = y * tot->GetBinContent(i);
        float stat = Fun->IntegralError(bini+(i-1)*binw*2,bini+i*binw*2) /binw/2;
        float eror = TMath::Sqrt( syst*syst + stat*stat);
	MisMax->SetBinContent(i,y + eror);
	MisMax->SetBinError(i,ey);
	MisMin->SetBinContent(i,y - eror);
	MisMin->SetBinError(i,ey);
	if ( y > 0. ) {
	  MisMax->SetBinError(i,ey * (1. + eror/y ));
	  MisMin->SetBinError(i,ey * (1. - eror/y ));
	}
	BinMin[i] = bini+(i-1)*binw*2;
	BinMax[i] = bini+i*binw*2;
	Leff[i] = Fun->Integral(BinMin[i],BinMax[i]) /binw/2;;
       }

       MisMax->Fit("Gun","rvee0");
       Gun->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LeffMax[i] = Gun->Integral(BinMin[i],BinMax[i]) / binw/2;
       }

       MisMin->Fit("Hun","rvee0");
       Hun->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LeffMin[i] = Hun->Integral(BinMin[i],BinMax[i]) / binw/2;
	LeffErr[i] = (LeffMax[i] - LeffMin[i]) / 2.;
       }

 //include text
  TPaveText *abc = new TPaveText(0.22,0.69,0.60,0.89,"brNDC");   
  abc->SetBorderSize(0);	
  abc->SetFillStyle(0);  
  abc->SetTextAlign(13); 
  abc->SetTextFont(42);	
  abc->SetTextSize(0.07);   
  TText *text = abc->AddText("SSVHPT tagger");   
  TText *text = abc->AddText("p_{T}(jet) > 30 GeV");   
  abc->Draw(""); 

  //include the official CMS label
  TPaveText *pt = new TPaveText(0.49,0.94,0.99,0.99,"brNDC");   
  pt->SetBorderSize(0);   
  pt->SetFillStyle(0);  
  pt->SetTextAlign(13);   
  pt->SetTextFont(42);   
  pt->SetTextSize(0.05);   
  TText *text = pt->AddText("CMS Preliminary 2010");   
  pt->Draw(""); 

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Data/MC Scale Factor

 TF1 *FUn = new TF1("FUn","([0]+[1]*x)",bini,bini+nbin*binw);
//  TF1 *FUn = new TF1("FUn","([0]+[1]*x)/(1+[2]*x)",bini,bini+nbin*binw);
//  TF1 *FUn = new TF1("FUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      FUn->SetLineColor(1);
      FUn->SetLineStyle(1);
      FUn->SetLineWidth(0.5);
 TF1 *GUn = new TF1("GUn","([0]+[1]*x)",bini,bini+nbin*binw);
//  TF1 *GUn = new TF1("GUn","([0]+[1]*x)/(1+[2]*x)",bini,bini+nbin*binw);
//  TF1 *GUn = new TF1("GUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      GUn->SetLineColor(1);
      GUn->SetLineStyle(2);
 TF1 *HUn = new TF1("HUn","([0]+[1]*x)",bini,bini+nbin*binw);
//  TF1 *HUn = new TF1("HUn","([0]+[1]*x)/(1+[2]*x)",bini,bini+nbin*binw);
//  TF1 *HUn = new TF1("HUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",bini,bini+nbin*binw);
      HUn->SetLineColor(1);
      HUn->SetLineStyle(2);

 TH1F* SFlight= new TH1F("SFlight","",nbin,bini,bini+nbin*binw);  // Scale Factor
       SFlight->Divide(Mistag,h4,1,1); 
 TH1F* SFneg= new TH1F("SFneg","",nbin,bini,bini+nbin*binw);  // neg tag Data/MC
       SFneg->Divide(Negtag,h5,1,1); 
       for (int i=1; i<nbin+1; i++) {
	float ey = SFneg->GetBinError(i);
	SFlight->SetBinError(i,ey);
       }

 TH1F* SF    = new TH1F("SF"  ,"",nbin/2,bini,bini+nbin*binw);
 TH1F* SFmax= new TH1F("SFmax","",nbin/2,bini,bini+nbin*binw);  // SF up
 TH1F* SFmin= new TH1F("SFmin","",nbin/2,bini,bini+nbin*binw);  // SF down
       for (int i=2; i<nbin+1; i = i+2) {
        float y  = SFlight->GetBinContent(i-1) + SFlight->GetBinContent(i);
        float ey = SFlight->GetBinError(i-1)*SFlight->GetBinError(i-1)
	         + SFlight->GetBinError(i)*SFlight->GetBinError(i);
	y = y / 2.;
	ey = TMath::Sqrt(ey) / 2.;
	int j = i/2;
        SF->SetBinContent(j,y);
        SF->SetBinError(j,ey);
       }

     pad4->cd();
       SF->Draw("E"); 
       SF->SetLineColor(2);
       SF->SetMarkerStyle(20);
       SF->SetMarkerColor(2);
       SF->SetMarkerSize(1.);
       SF->SetMinimum(0.0); 
       SF->SetMaximum(3.0);
       SF->GetXaxis()->SetNdivisions(509);
       SF->GetYaxis()->SetNdivisions(509);
       SF->SetTitleFont(42, "XYZ");
       SF->GetXaxis()->SetTitle(xtitle);
       SF->GetYaxis()->SetTitle("Data/MC");
       SF->GetYaxis()->SetTitleOffset(1.2);
       SF->SetTickLength(0.03, "XYZ");
       SF->SetLabelOffset(0.007,"XYZ");
       SF->SetLabelSize(0.06, "XYZ");
       SF->SetLabelFont(42, "XYZ");
       SF->SetTitleSize(0.06, "XYZ");
       SF->SetNdivisions(509,"XYZ");
       SF->Draw("Esame"); 
       SF->Fit("FUn","rvee"); 
       FUn->Draw("same"); 
 
       for (int i=1; i<nbin/2+1; i++) {
	float y = SF->GetBinContent(i);
	float ey = SF->GetBinError(i);
	float syst = y * tot->GetBinContent(i);
        float eror = TMath::Sqrt( syst*syst + ey*ey);
        float yfit = FUn->Integral(bini+(i-1)*binw*2,bini+i*binw*2) / binw/2;
	SFmax->SetBinContent(i,yfit + eror);
	SFmax->SetBinError(i,ey);
	SFmin->SetBinContent(i,yfit - eror);
	SFmin->SetBinError(i,ey);
	Lsf[i] = FUn->Integral(BinMin[i],BinMax[i]) / binw/2;
       }

       SFmax->Fit("GUn","rvee0");
       GUn->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LsfMax[i] = GUn->Integral(BinMin[i],BinMax[i]) / binw/2;
       }

       SFmin->Fit("HUn","rvee0");
       HUn->Draw("same"); 
       for (int i=1; i<nbin/2+1; i++) {
	LsfMin[i] = HUn->Integral(BinMin[i],BinMax[i]) / binw/2;
	LsfErr[i] = (LsfMax[i] - LsfMin[i]) / 2.;
       }

  //include text
  TPaveText *abc = new TPaveText(0.22,0.69,0.60,0.89,"brNDC");   
  abc->SetBorderSize(0);	
  abc->SetFillStyle(0);  
  abc->SetTextAlign(13); 
  abc->SetTextFont(42);	
  abc->SetTextSize(0.07);   
  TText *text = abc->AddText("SSVHPT tagger");   
  TText *text = abc->AddText("p_{T}(jet) > 30 GeV");   
  abc->Draw(""); 

  //include the official CMS label
  TPaveText *pt = new TPaveText(0.49,0.94,0.99,0.99,"brNDC");   
  pt->SetBorderSize(0);   
  pt->SetFillStyle(0);  
  pt->SetTextAlign(13);   
  pt->SetTextFont(42);   
  pt->SetTextSize(0.05);   
  TText *text = pt->AddText("CMS Preliminary 2010");   
  pt->Draw(""); 

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  c1->Update();
}

