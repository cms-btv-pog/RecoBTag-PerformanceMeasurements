#include <iostream>
#include <TROOT.h>

plot()
{
// *****************************************************************************

 Int_t stati=0;
 Bool_t  fit=0;
 Bool_t logy=0;

 float sfK0 = 1.95;
 float sfLA = 2.2;

 char* hname  = "hData_Tagger"; 
 char* hlight = "hLightFlav_Tagger"; char* hcflav = "hCFlav_Tagger"; char* hnoflav = "hNoFlav_Tagger"; 

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
pad1->SetLogy(1);
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.05);
   pad1->SetLeftMargin(0.15);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(1);
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
pad4->SetLogy(0);
   pad4->SetTopMargin(0.05);
   pad4->SetBottomMargin(0.15);
   pad4->SetRightMargin(0.05);
   pad4->SetLeftMargin(0.15);

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

    pad1->cd();
 TFile *f1 = new TFile("../output/Data_ptGT30_TCHEL.root");
 TFile *g1 = new TFile("../output/MC_QCD-Pt15_ptGT30_TCHEL.root");
 char* htitle = "TCHE discriminator"; 
 int nbin = 100; float binmin = -25, binmax = 25;

 TH1F* h0= new TH1F("h0","", nbin, binmin, binmax);
 TH1F* h1= new TH1F("h1","", nbin, binmin, binmax);
 TH1F* h2= new TH1F("h2","", nbin, binmin, binmax);
 TH1F* hl= new TH1F("hl","", nbin, binmin, binmax);
 TH1F* hc= new TH1F("hc","", nbin, binmin, binmax);
 TH1F* hno= new TH1F("hno","", nbin, binmin, binmax);
 TH1F* hv0= new TH1F("hv0","", nbin, binmin, binmax);
       f1->cd();
 TH1F* h10 = (TH1F*)gROOT->FindObject(hname);
       h1->Add(h10,h1,1,0);
       g1->cd();
 TH1F* h20 = (TH1F*)gROOT->FindObject(hname);
 TH1F* hks = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_K0s");
 TH1F* hla = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_Lam");
       h2->Add(h20,hks,1,sfK0-1);
       h2->Add(h2,hla,1,sfLA-1);
       float norm = h1->Integral(0,nbin+1)/h2->Integral(0,nbin+1);
       h2->Add(h2,h2,norm,0);

 TH1F* hl0 = (TH1F*)gROOT->FindObject(hlight);
       h0->Add(hl0,h0,norm,0);
       hl->Add(hks,hla,sfK0-1,sfLA-1);
       hl->Add(hl0,hl,norm,norm);
 TH1F* hno0 = (TH1F*)gROOT->FindObject(hnoflav);
       hno->Add(hno0,hl,norm,1);
 TH1F* hc0 = (TH1F*)gROOT->FindObject(hcflav);
       hc->Add(hc0,hno,norm,1);
       hv0->Add(hks,hla,sfK0*norm,sfLA*norm);

 TH1F* H1= new TH1F("H1","", nbin, binmin, binmax);
 TH1F* H2= new TH1F("H2","", nbin, binmin, binmax);
 TH1F* H3= new TH1F("H3","", nbin, binmin, binmax);
       for (int i=1; i<nbin/2+1; i++) {
         H1->SetBinContent(i,h2->GetBinContent(i));  
         H2->SetBinContent(i,hc->GetBinContent(i));  
         H3->SetBinContent(i,hno->GetBinContent(i));
       }  

       h1->Draw("E"); 
       h1->SetMarkerStyle(20);h1->SetMarkerSize(0.5); 
       h1->SetMarkerColor(kBlack);
       h1->SetLineColor(kBlack);
       h1->SetLineStyle(1);
       h1->SetLineWidth(1);
//        hv->Draw("Hsame"); 
//        hv->SetFillColor(kOrange);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kRed-4);
       H1->Draw("Hsame"); 
       H1->SetFillColor(kRed);
       hc->Draw("Hsame"); 
       hc->SetFillColor(kGreen-4);
       H2->Draw("Hsame"); 
       H2->SetFillColor(kGreen);
       hno->Draw("Hsame"); 
       hno->SetFillColor(kBlue-4);
       H3->Draw("Hsame"); 
       H3->SetFillColor(kBlue);
//        hl->Draw("Hsame"); 
//        hl->SetFillColor(kBlue-7);
//        hv0->Draw("Hsame"); 
//        hv0->SetFillColor(kBlue-7);
//        h0->Draw("Hsame"); 
//        h0->SetFillColor(kBlue);
       h1->SetTickLength(0.03, "YZ");
       h1->SetTickLength(-0.03,"X");
       h1->SetLabelOffset(0.023,"X");
       h1->SetLabelOffset(0.007,"Y");
       h1->SetLabelSize(0.07, "XYZ");
       h1->SetLabelFont(42, "XYZ");
       h1->SetTitleSize(0.07, "XYZ");
       h1->GetXaxis()->SetTitle(htitle);
       h1->GetYaxis()->SetTitle("mistag rate");
       h1->SetNdivisions(509,"XYZ");
       h1->Draw("Esame"); 
//        h1->SetMinimum(1); 
//        h1->SetMaximum(3000000);  

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

    pad2->cd();
 TFile *f1 = new TFile("../output/Data_ptGT30_TCHPT.root");
 TFile *g1 = new TFile("../output/MC_QCD-Pt15_ptGT30_TCHPT.root");
 char* htitle = "TCHP discriminator"; 
 int nbin = 100; float binmin = -25, binmax = 25;

 TH1F* h0= new TH1F("h0","", nbin, binmin, binmax);
 TH1F* h1= new TH1F("h1","", nbin, binmin, binmax);
 TH1F* h2= new TH1F("h2","", nbin, binmin, binmax);
 TH1F* hl= new TH1F("hl","", nbin, binmin, binmax);
 TH1F* hc= new TH1F("hc","", nbin, binmin, binmax);
 TH1F* hno= new TH1F("hno","", nbin, binmin, binmax);
 TH1F* hv0= new TH1F("hv0","", nbin, binmin, binmax);
       f1->cd();
 TH1F* h10 = (TH1F*)gROOT->FindObject(hname);
       h1->Add(h10,h1,1,0);
       g1->cd();
 TH1F* h20 = (TH1F*)gROOT->FindObject(hname);
 TH1F* hks = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_K0s");
 TH1F* hla = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_Lam");
       h2->Add(h20,hks,1,sfK0-1);
       h2->Add(h2,hla,1,sfLA-1);
       float norm = h1->Integral(0,nbin+1)/h2->Integral(0,nbin+1);
       h2->Add(h2,h2,norm,0);

 TH1F* hl0 = (TH1F*)gROOT->FindObject(hlight);
       h0->Add(hl0,h0,norm,0);
       hl->Add(hks,hla,sfK0-1,sfLA-1);
       hl->Add(hl0,hl,norm,norm);
 TH1F* hno0 = (TH1F*)gROOT->FindObject(hnoflav);
       hno->Add(hno0,hl,norm,1);
 TH1F* hc0 = (TH1F*)gROOT->FindObject(hcflav);
       hc->Add(hc0,hno,norm,1);
       hv0->Add(hks,hla,sfK0*norm,sfLA*norm);

 TH1F* H1= new TH1F("H1","", nbin, binmin, binmax);
 TH1F* H2= new TH1F("H2","", nbin, binmin, binmax);
 TH1F* H3= new TH1F("H3","", nbin, binmin, binmax);
       for (int i=1; i<nbin/2+1; i++) {
         H1->SetBinContent(i,h2->GetBinContent(i));  
         H2->SetBinContent(i,hc->GetBinContent(i));  
         H3->SetBinContent(i,hno->GetBinContent(i));
       }  

       h1->Draw("E"); 
       h1->SetMarkerStyle(20);h1->SetMarkerSize(0.5); 
       h1->SetMarkerColor(kBlack);
       h1->SetLineColor(kBlack);
       h1->SetLineStyle(1);
       h1->SetLineWidth(1);
//        hv->Draw("Hsame"); 
//        hv->SetFillColor(kOrange);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kRed-4);
       H1->Draw("Hsame"); 
       H1->SetFillColor(kRed);
       hc->Draw("Hsame"); 
       hc->SetFillColor(kGreen-4);
       H2->Draw("Hsame"); 
       H2->SetFillColor(kGreen);
       hno->Draw("Hsame"); 
       hno->SetFillColor(kBlue-4);
       H3->Draw("Hsame"); 
       H3->SetFillColor(kBlue);
//        hl->Draw("Hsame"); 
//        hl->SetFillColor(kBlue-7);
//        hv0->Draw("Hsame"); 
//        hv0->SetFillColor(kBlue-7);
//        h0->Draw("Hsame"); 
//        h0->SetFillColor(kBlue);
       h1->SetTickLength(0.03, "YZ");
       h1->SetTickLength(-0.03,"X");
       h1->SetLabelOffset(0.023,"X");
       h1->SetLabelOffset(0.007,"Y");
       h1->SetLabelSize(0.07, "XYZ");
       h1->SetLabelFont(42, "XYZ");
       h1->SetTitleSize(0.07, "XYZ");
       h1->GetXaxis()->SetTitle(htitle);
       h1->SetNdivisions(509,"XYZ");
       h1->Draw("Esame"); 
//        h1->SetMinimum(1); 
//        h1->SetMaximum(3000000);  

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

    pad3->cd();
 TFile *f1 = new TFile("../output/Data_ptGT30_SSVHEM.root");
 TFile *g1 = new TFile("../output/MC_QCD-Pt15_ptGT30_SSVHEM.root");
 char* htitle = "SSVHE discriminator"; 
 int nbin = 100; float binmin = -5, binmax = 5;

 TH1F* h0= new TH1F("h0","", nbin, binmin, binmax);
 TH1F* h1= new TH1F("h1","", nbin, binmin, binmax);
 TH1F* h2= new TH1F("h2","", nbin, binmin, binmax);
 TH1F* hl= new TH1F("hl","", nbin, binmin, binmax);
 TH1F* hc= new TH1F("hc","", nbin, binmin, binmax);
 TH1F* hno= new TH1F("hno","", nbin, binmin, binmax);
 TH1F* hv0= new TH1F("hv0","", nbin, binmin, binmax);
       f1->cd();
 TH1F* h10 = (TH1F*)gROOT->FindObject(hname);
       h1->Add(h10,h1,1,0);
       g1->cd();
 TH1F* h20 = (TH1F*)gROOT->FindObject(hname);
 TH1F* hks = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_K0s");
 TH1F* hla = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_Lam");
       h2->Add(h20,hks,1,sfK0-1);
       h2->Add(h2,hla,1,sfLA-1);
       float norm = h1->Integral(0,nbin+1)/h2->Integral(0,nbin+1);
       h2->Add(h2,h2,norm,0);

 TH1F* hl0 = (TH1F*)gROOT->FindObject(hlight);
       h0->Add(hl0,h0,norm,0);
       hl->Add(hks,hla,sfK0-1,sfLA-1);
       hl->Add(hl0,hl,norm,norm);
 TH1F* hno0 = (TH1F*)gROOT->FindObject(hnoflav);
       hno->Add(hno0,hl,norm,1);
 TH1F* hc0 = (TH1F*)gROOT->FindObject(hcflav);
       hc->Add(hc0,hno,norm,1);
       hv0->Add(hks,hla,sfK0*norm,sfLA*norm);

 TH1F* H1= new TH1F("H1","", nbin, binmin, binmax);
 TH1F* H2= new TH1F("H2","", nbin, binmin, binmax);
 TH1F* H3= new TH1F("H3","", nbin, binmin, binmax);
       for (int i=1; i<nbin/2+1; i++) {
         H1->SetBinContent(i,h2->GetBinContent(i));  
         H2->SetBinContent(i,hc->GetBinContent(i));  
         H3->SetBinContent(i,hno->GetBinContent(i));
       }  

       h1->Draw("E"); 
       h1->SetMarkerStyle(20);h1->SetMarkerSize(0.5); 
       h1->SetMarkerColor(kBlack);
       h1->SetLineColor(kBlack);
       h1->SetLineStyle(1);
       h1->SetLineWidth(1);
//        hv->Draw("Hsame"); 
//        hv->SetFillColor(kOrange);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kRed-4);
       H1->Draw("Hsame"); 
       H1->SetFillColor(kRed);
       hc->Draw("Hsame"); 
       hc->SetFillColor(kGreen-4);
       H2->Draw("Hsame"); 
       H2->SetFillColor(kGreen);
       hno->Draw("Hsame"); 
       hno->SetFillColor(kBlue-4);
       H3->Draw("Hsame"); 
       H3->SetFillColor(kBlue);
//        hl->Draw("Hsame"); 
//        hl->SetFillColor(kBlue-7);
//        hv0->Draw("Hsame"); 
//        hv0->SetFillColor(kBlue-7);
//        h0->Draw("Hsame"); 
//        h0->SetFillColor(kBlue);
       h1->SetTickLength(0.03, "YZ");
       h1->SetTickLength(-0.03,"X");
       h1->SetLabelOffset(0.023,"X");
       h1->SetLabelOffset(0.007,"Y");
       h1->SetLabelSize(0.07, "XYZ");
       h1->SetLabelFont(42, "XYZ");
       h1->SetTitleSize(0.07, "XYZ");
       h1->GetXaxis()->SetTitle(htitle);
       h1->SetNdivisions(509,"XYZ");
       h1->Draw("Esame"); 
//        h1->SetMinimum(1); 
//        h1->SetMaximum(3000000);  

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

    pad4->cd();
 TFile *f1 = new TFile("../output/Data_ptGT30_SSVHPT.root");
 TFile *g1 = new TFile("../output/MC_QCD-Pt15_ptGT30_SSVHPT.root");
 char* htitle = "SSVHP discriminator"; 
 int nbin = 100; float binmin = -5, binmax = 5;

 TH1F* h0= new TH1F("h0","", nbin, binmin, binmax);
 TH1F* h1= new TH1F("h1","", nbin, binmin, binmax);
 TH1F* h2= new TH1F("h2","", nbin, binmin, binmax);
 TH1F* hl= new TH1F("hl","", nbin, binmin, binmax);
 TH1F* hc= new TH1F("hc","", nbin, binmin, binmax);
 TH1F* hno= new TH1F("hno","", nbin, binmin, binmax);
 TH1F* hv0= new TH1F("hv0","", nbin, binmin, binmax);
       f1->cd();
 TH1F* h10 = (TH1F*)gROOT->FindObject(hname);
       h1->Add(h10,h1,1,0);
       g1->cd();
 TH1F* h20 = (TH1F*)gROOT->FindObject(hname);
 TH1F* hks = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_K0s");
 TH1F* hla = (TH1F*)gROOT->FindObject("hAllFlav_Tagger_Lam");
       h2->Add(h20,hks,1,sfK0-1);
       h2->Add(h2,hla,1,sfLA-1);
       float norm = h1->Integral(0,nbin+1)/h2->Integral(0,nbin+1);
       h2->Add(h2,h2,norm,0);

 TH1F* hl0 = (TH1F*)gROOT->FindObject(hlight);
       h0->Add(hl0,h0,norm,0);
       hl->Add(hks,hla,sfK0-1,sfLA-1);
       hl->Add(hl0,hl,norm,norm);
 TH1F* hno0 = (TH1F*)gROOT->FindObject(hnoflav);
       hno->Add(hno0,hl,norm,1);
 TH1F* hc0 = (TH1F*)gROOT->FindObject(hcflav);
       hc->Add(hc0,hno,norm,1);
       hv0->Add(hks,hla,sfK0*norm,sfLA*norm);

 TH1F* H1= new TH1F("H1","", nbin, binmin, binmax);
 TH1F* H2= new TH1F("H2","", nbin, binmin, binmax);
 TH1F* H3= new TH1F("H3","", nbin, binmin, binmax);
       for (int i=1; i<nbin/2+1; i++) {
         H1->SetBinContent(i,h2->GetBinContent(i));  
         H2->SetBinContent(i,hc->GetBinContent(i));  
         H3->SetBinContent(i,hno->GetBinContent(i));
       }  

       h1->Draw("E"); 
       h1->SetMarkerStyle(20);h1->SetMarkerSize(0.5); 
       h1->SetMarkerColor(kBlack);
       h1->SetLineColor(kBlack);
       h1->SetLineStyle(1);
       h1->SetLineWidth(1);
//        hv->Draw("Hsame"); 
//        hv->SetFillColor(kOrange);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kRed-4);
       H1->Draw("Hsame"); 
       H1->SetFillColor(kRed);
       hc->Draw("Hsame"); 
       hc->SetFillColor(kGreen-4);
       H2->Draw("Hsame"); 
       H2->SetFillColor(kGreen);
       hno->Draw("Hsame"); 
       hno->SetFillColor(kBlue-4);
       H3->Draw("Hsame"); 
       H3->SetFillColor(kBlue);
//        hl->Draw("Hsame"); 
//        hl->SetFillColor(kBlue-7);
//        hv0->Draw("Hsame"); 
//        hv0->SetFillColor(kBlue-7);
//        h0->Draw("Hsame"); 
//        h0->SetFillColor(kBlue);
       h1->SetTickLength(0.03, "YZ");
       h1->SetTickLength(-0.03,"X");
       h1->SetLabelOffset(0.023,"X");
       h1->SetLabelOffset(0.007,"Y");
       h1->SetLabelSize(0.07, "XYZ");
       h1->SetLabelFont(42, "XYZ");
       h1->SetTitleSize(0.07, "XYZ");
       h1->GetXaxis()->SetTitle(htitle);
       h1->SetNdivisions(509,"XYZ");
       h1->Draw("Esame"); 
//        h1->SetMinimum(1); 
//        h1->SetMaximum(3000000);  

  TLegend* leg = new TLegend(0.20,0.63,0.50,0.93);
//     leg->SetHeader("jet pt > 30 GeV");
//     leg->SetMargin(0.12);
//     leg->SetTextSize(0.035);
//     leg->SetFillStyle(1);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42);
    leg->AddEntry(h1," Data","PL");
    leg->AddEntry(h2," MC b","F");
    leg->AddEntry(hc," MC c","F");
//     leg->AddEntry(hl," V^{0} reweighting","F");
//     leg->AddEntry(h0," MC light","F");
    leg->AddEntry(hno," MC udsg","F");
    leg->Draw();
    
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

