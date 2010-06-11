#include <iostream>
#include <TROOT.h>
#include <TPaveText.h>

plot()
{
 TFile *f1 = new TFile("../output/MC_QCD-Pt15_ptGT30_SSVHET.root");

// *****************************************************************************

 Int_t stati=0;
 Bool_t  fit=0;
 Bool_t logy=0;

// *****************************************************************************

 TCanvas *c1 = new TCanvas("c1", "plots",200,10,700,730);
 c1->SetFillColor(10);
 c1->SetFillStyle(4000);
 c1->SetBorderSize(2);

// *****************************************************************************

// TPaveLabel *p01 = new TPaveLabel(0.05,0.93,0.95,0.97,
// //                   "SSVHE in minbias udsg jets (MC CMSSW 3_5_7)","br");
//                   "SSVHE in udsg jets (/MinBias/Spring10-START3X_V26A_356ReReco-v1/)","br");
// //                   "SSVHE in udsg jets (/QCD_Pt30/Spring10-START3X_V26_S09-v1)","br");
// // p01->SetFillColor(7);
// p01->SetFillColor(0);
// p01->SetFillStyle(3017);
// p01->SetTextSize(0.8);
// p01->Draw();

pad1 = new TPad("pad1","This is pad1",0.03,0.49,0.33,0.91,21);
pad2 = new TPad("pad2","This is pad2",0.36,0.49,0.66,0.91,21);
pad3 = new TPad("pad3","This is pad2",0.69,0.49,0.99,0.91,21);
pad4 = new TPad("pad4","This is pad4",0.03,0.04,0.33,0.46,21);
pad5 = new TPad("pad5","This is pad5",0.36,0.04,0.66,0.46,21);
pad6 = new TPad("pad6","This is pad6",0.69,0.04,0.99,0.46,21);

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
pad2->SetLogy(logy);
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.15);
   pad2->SetRightMargin(0.05);
   pad2->SetLeftMargin(0.15);

pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(logy);
   pad3->SetTopMargin(0.05);
   pad3->SetBottomMargin(0.15);
   pad3->SetRightMargin(0.05);
   pad3->SetLeftMargin(0.15);

pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(logy);
   pad4->SetTopMargin(0.05);
   pad4->SetBottomMargin(0.15);
   pad4->SetRightMargin(0.05);
   pad4->SetLeftMargin(0.15);

pad5->SetFillColor(0);
pad5->SetBorderMode(0);
pad5->SetFrameFillColor(10);
pad5->Draw();
pad5->SetLogy(logy);
   pad5->SetTopMargin(0.05);
   pad5->SetBottomMargin(0.15);
   pad5->SetRightMargin(0.05);
   pad5->SetLeftMargin(0.15);

pad6->SetFillColor(0);
pad6->SetBorderMode(0);
pad6->SetFrameFillColor(10);
pad6->Draw();
pad6->SetLogy(logy);
   pad6->SetTopMargin(0.05);
   pad6->SetBottomMargin(0.15);
   pad6->SetRightMargin(0.05);
   pad6->SetLeftMargin(0.15);

//$$ gStyle->SetOptDate(1);
gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleColor(29);
gStyle->SetTitleW(0.9);
gStyle->SetTitleH(0.1);
gStyle->SetOptStat(stati);

if (fit) {
gStyle->SetOptFit(111);
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.1);
} else {
gStyle->SetOptFit(0);
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.1);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 f1->cd();

 TH1F* g00 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* g01 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Gam");
 TH1F* g02 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_K0s");
 TH1F* g03 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Lam");
 TH1F* g04 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Bwd");
 TH1F* g05 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Cwd");
 TH1F* g06 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Tau");
 TH1F* g07 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Int");
 TH1F* g08 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Fak");
 TH1F* g09 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Bad");
 TH1F* g10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Oth");

 pad1->cd();

 TH1F* h0 = new TH1F("h0","",100,-5.,5.);
 TH1F* h1 = new TH1F("h1","",50,0.,5.);
 TH1F* h2 = new TH1F("h2","",50,0.,5.);
       h0->Add(g01,h0,1,0);

       for (Int_t i=1; i<51; i++) {
         h1->SetBinContent(i,h0->GetBinContent(i+50));  
         h2->SetBinContent(i,h0->GetBinContent(51-i));  
       }  
       h1->Draw("E"); 
       h1->SetLineColor(2);
       h1->SetLineStyle(1);
       h1->SetLineWidth(2);
       h1->GetXaxis()->SetLabelSize(0.08);
       h1->GetYaxis()->SetLabelSize(0.08);
       h1->GetXaxis()->SetTitleSize(0.08);
       h1->GetXaxis()->SetTitle("Conversion");
       h1->GetXaxis()->SetTitleColor(1);
       h1->GetXaxis()->SetNdivisions(509);
       h1->GetYaxis()->SetNdivisions(509);
       h2->Draw("same"); h1->Draw("Esame"); 
       h2->SetLineColor(4);
       h2->SetLineStyle(1);
       h2->SetLineWidth(2);

  TLegend* leg = new TLegend(0.48,0.64,0.90,0.92);
    leg->SetHeader("SSVHE");
    leg->AddEntry(h1," pos tag","L");
    leg->AddEntry(h2," neg tag","L");
    leg->Draw();

 pad2->cd();

 TH1F* h0 = new TH1F("h0","",100,-5.,5.);
 TH1F* h1 = new TH1F("h1","",50,0.,5.);
 TH1F* h2 = new TH1F("h2","",50,0.,5.);
       h0->Add(g02,h0,1,0);

       for (Int_t i=1; i<51; i++) {
         h1->SetBinContent(i,h0->GetBinContent(i+50));  
         h2->SetBinContent(i,h0->GetBinContent(51-i));  
       }  
       h1->Draw("E"); 
       h1->SetLineColor(2);
       h1->SetLineStyle(1);
       h1->SetLineWidth(2);
       h1->GetXaxis()->SetLabelSize(0.08);
       h1->GetYaxis()->SetLabelSize(0.08);
       h1->GetXaxis()->SetTitleSize(0.08);
       h1->GetXaxis()->SetTitle("K0s");
       h1->GetXaxis()->SetTitleColor(1);
       h1->GetXaxis()->SetNdivisions(509);
       h1->GetYaxis()->SetNdivisions(509);
       h2->Draw("same"); h1->Draw("Esame"); 
       h2->SetLineColor(4);
       h2->SetLineStyle(1);
       h2->SetLineWidth(2);

 pad3->cd();

 TH1F* h0 = new TH1F("h0","",100,-5.,5.);
 TH1F* h1 = new TH1F("h1","",50,0.,5.);
 TH1F* h2 = new TH1F("h2","",50,0.,5.);
       h0->Add(g03,h0,1,0);

       for (Int_t i=1; i<51; i++) {
         h1->SetBinContent(i,h0->GetBinContent(i+50));  
         h2->SetBinContent(i,h0->GetBinContent(51-i));  
       }  
       h1->Draw("E"); 
       h1->SetLineColor(2);
       h1->SetLineStyle(1);
       h1->SetLineWidth(2);
       h1->GetXaxis()->SetLabelSize(0.08);
       h1->GetYaxis()->SetLabelSize(0.08);
       h1->GetXaxis()->SetTitleSize(0.08);
       h1->GetXaxis()->SetTitle("#Lambda");
       h1->GetXaxis()->SetTitleColor(1);
       h1->GetXaxis()->SetNdivisions(509);
       h1->GetYaxis()->SetNdivisions(509);
       h2->Draw("same"); h1->Draw("Esame"); 
       h2->SetLineColor(4);
       h2->SetLineStyle(1);
       h2->SetLineWidth(2);

 pad4->cd();

 TH1F* h0 = new TH1F("h0","",100,-5.,5.);
 TH1F* h1 = new TH1F("h1","",50,0.,5.);
 TH1F* h2 = new TH1F("h2","",50,0.,5.);
       h0->Add(g07,h0,1,0);

       for (Int_t i=1; i<51; i++) {
         h1->SetBinContent(i,h0->GetBinContent(i+50));  
         h2->SetBinContent(i,h0->GetBinContent(51-i));  
       }  
       h1->Draw("E"); 
       h1->SetLineColor(2);
       h1->SetLineStyle(1);
       h1->SetLineWidth(2);
       h1->GetXaxis()->SetLabelSize(0.08);
       h1->GetYaxis()->SetLabelSize(0.08);
       h1->GetXaxis()->SetTitleSize(0.08);
//        h1->GetXaxis()->SetTitle("Discriminator x10-10");
       h1->GetXaxis()->SetTitle("hadronic int.");
       h1->GetXaxis()->SetTitleColor(1);
       h1->GetXaxis()->SetNdivisions(509);
       h1->GetYaxis()->SetNdivisions(509);
//        h1->SetMinimum(20); h1->SetMaximum(3000000);  
       h2->Draw("same"); h1->Draw("Esame"); 
       h2->SetLineColor(4);
       h2->SetLineStyle(1);
       h2->SetLineWidth(2);

 pad5->cd();

 TH1F* h0 = new TH1F("h0","",100,-5.,5.);
 TH1F* h1 = new TH1F("h1","",50,0.,5.);
 TH1F* h2 = new TH1F("h2","",50,0.,5.);
       h0->Add(g08,h0,1,0);

       for (Int_t i=1; i<51; i++) {
         h1->SetBinContent(i,h0->GetBinContent(i+50));  
         h2->SetBinContent(i,h0->GetBinContent(51-i));  
       }  
       h1->Draw("E"); 
       h1->SetLineColor(2);
       h1->SetLineStyle(1);
       h1->SetLineWidth(2);
       h1->GetXaxis()->SetLabelSize(0.08);
       h1->GetYaxis()->SetLabelSize(0.08);
       h1->GetXaxis()->SetTitleSize(0.08);
//        h1->GetXaxis()->SetTitle("Discriminator x10-10");
       h1->GetXaxis()->SetTitle("Fake");
       h1->GetXaxis()->SetTitleColor(1);
       h1->GetXaxis()->SetNdivisions(509);
       h1->GetYaxis()->SetNdivisions(509);
//        h1->SetMinimum(20); h1->SetMaximum(3000000);  
       h2->Draw("same"); h1->Draw("Esame"); 
       h2->SetLineColor(4);
       h2->SetLineStyle(1);
       h2->SetLineWidth(2);

 pad6->cd();

 TH1F* h0 = new TH1F("h0","",100,-5.,5.);
 TH1F* h1 = new TH1F("h1","",50,0.,5.);
 TH1F* h2 = new TH1F("h2","",50,0.,5.);
       h0->Add(g09,g10,1,1);

       for (Int_t i=1; i<51; i++) {
         h1->SetBinContent(i,h0->GetBinContent(i+50));  
         h2->SetBinContent(i,h0->GetBinContent(51-i));  
       }  
       h1->Draw("E"); 
       h1->SetLineColor(2);
       h1->SetLineStyle(1);
       h1->SetLineWidth(2);
       h1->GetXaxis()->SetLabelSize(0.08);
       h1->GetYaxis()->SetLabelSize(0.08);
       h1->GetXaxis()->SetTitleSize(0.08);
//        h1->GetXaxis()->SetTitle("Discriminator x10-10");
       h1->GetXaxis()->SetTitle("Others");
       h1->GetXaxis()->SetTitleColor(1);
       h1->GetXaxis()->SetNdivisions(509);
       h1->GetYaxis()->SetNdivisions(509);
//        h1->SetMinimum(20); h1->SetMaximum(3000000);  
       h2->Draw("same"); h1->Draw("Esame"); 
       h2->SetLineColor(4);
       h2->SetLineStyle(1);
       h2->SetLineWidth(2);

float n1 = hLightFlav_PosTag_Gam->GetEntries();
float n2 = hLightFlav_PosTag_K0Lam->GetEntries();
float n4 = hLightFlav_PosTag_BCTau->GetEntries();
float n7 = hLightFlav_PosTag_Int->GetEntries();
float n8 = hLightFlav_PosTag_Fak->GetEntries();
float n10= hLightFlav_PosTag_Bad->GetEntries() - hLightFlav_NegTag_Bad->GetEntries();
    n10 += hLightFlav_PosTag_Oth->GetEntries() - hLightFlav_NegTag_Oth->GetEntries();
float n0 = hLightFlav_PosTag_JetPt->GetEntries() -n1 -n2 -n4 -n7 -n8 -n10;

cout << " Light mistags : " << hLightFlav_PosTag_JetPt->GetEntries() / hData_JetPt->GetEntries() << endl;
cout << " Conversion :  " << n1 / hLightFlav_PosTag_JetPt->GetEntries() << endl;
cout << " K0s+Lambda :  " << n2 / hLightFlav_PosTag_JetPt->GetEntries() << endl;
cout << " B+C+tau :     " << n4 / hLightFlav_PosTag_JetPt->GetEntries() << endl;
cout << " 2ndary int.:  " << n7 / hLightFlav_PosTag_JetPt->GetEntries() << endl;
cout << " Fake :        " << n8 / hLightFlav_PosTag_JetPt->GetEntries() << endl;
cout << " Oth Pos-Neg : " << n10/ hLightFlav_PosTag_JetPt->GetEntries() << endl;
cout << " Others :      " << n0 / hLightFlav_PosTag_JetPt->GetEntries() << endl;

 float SFbad = 1 + n10 / n8;
cout << endl;
cout << " Fakes =          " << hLightFlav_PosTag_Fak->GetEntries() << endl;
cout << " Others Pos-Neg = " << hLightFlav_PosTag_Oth->GetEntries()+hLightFlav_PosTag_Bad->GetEntries() 
                    << " - " << hLightFlav_NegTag_Oth->GetEntries()+hLightFlav_NegTag_Bad->GetEntries() << endl;
cout << " SFbad " << SFbad << endl;
cout << endl;

  c1->Update();
}

