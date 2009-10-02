#include <iostream>
#include <TROOT.h>
#include <TPaveText.h>

plot()
{
 TFile *f1 = new TFile("../output/QCDpt300_preprod_JPm.root");

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

pad1 = new TPad("pad1","This is pad1",0.01,0.01,0.99,0.99,21);

pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
if (logy) pad1->SetLogy();

//$$ gStyle->SetOptDate(1);
gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleColor(29);
gStyle->SetTitleW(0.5);
gStyle->SetTitleH(0.05);
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

    pad1->cd();

      f1->cd();
 TH1F* h0= new TH1F("h0","Jet Probability : QCD MC",100,-50.,50.);
 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* g1= new TH1F("g1","",100,-50.,50.);
 TH1F* g2= new TH1F("g2","",100,-50.,50.);
 TH1F* g3= new TH1F("g3","",100,-50.,50.);
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Veto_Tagger");
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Veto_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Veto_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Veto_Tagger");
       h0->Add(h00,h0,1,0);
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h2,1,1);
       for (int i=1; i<51; i++) {
         g1->SetBinContent(i,h1->GetBinContent(i));  
         g2->SetBinContent(i,h2->GetBinContent(i));  
         g3->SetBinContent(i,h3->GetBinContent(i));
       }  

       h0->Draw("E"); 
       h0->SetLineColor(1);
       h0->SetLineStyle(1);
       h0->SetLineWidth(1);
       h3->Draw("Hsame"); 
       h3->SetFillColor(kRed);
       g3->Draw("Hsame"); 
       g3->SetFillColor(kRed+1);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kBlue);
       g2->Draw("Hsame"); 
       g2->SetFillColor(kBlue+1);
       h1->Draw("Hsame"); 
       h1->SetFillColor(kGreen);
       g1->Draw("Hsame"); 
       g1->SetFillColor(kGreen+1);
       h0->SetTickLength(-0.02,"X");
       h0->SetLabelOffset(0.015,"X");
       h0->GetXaxis()->SetLabelSize(0.04);
       h0->GetYaxis()->SetLabelSize(0.04);
       h0->GetXaxis()->SetTitleSize(0.04);
       h0->GetXaxis()->SetTitle("Discriminator #times 20");
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h1->SetMaximum(3000000);  

  TLegend* leg = new TLegend(0.15,0.71,0.32,0.88);
//     leg->SetHeader(" all jets");
    leg->AddEntry(h3," b-jet","F");
    leg->AddEntry(h2," c-jet","F");
    leg->AddEntry(h1," udsg ","F");
    leg->Draw();

//   TPaveText* pt1 = new TPaveText(0.70,0.82,0.88,0.88,"NDC");
//   TText* t1 = pt1->AddText("QCD MC");
// //   TPaveText* pt1 = new TPaveText(0.55,0.82,0.88,0.88,"NDC");
// //   TText* t1 = pt1->AddText("QCD muon-jet MC");
//   t1->SetTextColor(1);
//   t1->SetTextSize(0.04);
//   pt1->Draw();

 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* k1= new TH1F("k1","",100,-50.,50.);
 TH1F* k2= new TH1F("k2","",100,-50.,50.);
 TH1F* k3= new TH1F("k3","",100,-50.,50.);
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h1,1,1);
       for (int i=1; i<51; i++) {
         k1->SetBinContent(i,h1->GetBinContent(i));
         k2->SetBinContent(i,h2->GetBinContent(i));  
         k3->SetBinContent(i,h3->GetBinContent(i));  
       }  
       k3->Draw("Hsame"); 
       k3->SetLineColor(kRed);
       k3->SetLineWidth(2);
       k2->Draw("Hsame"); 
       k2->SetLineColor(kBlue);
       k2->SetLineWidth(2);
       k1->Draw("Hsame"); 
       k1->SetLineColor(kGreen);
       k1->SetLineWidth(2);

  c1->Update();
}

