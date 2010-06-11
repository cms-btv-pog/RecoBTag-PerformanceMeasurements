#include <iostream>
#include <TROOT.h>

plot()
{

// *****************************************************************************

 int stati=0;
 bool  fit=0;
 bool logy=1;

 int nbin = 50;
 float SFk0 = 1.95 - 1.;
 float SFla = 2.2 - 1.;

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
pad3->SetLogy(0);
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

//$$ gStyle->SetOptDate(1);
gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleColor(29);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
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
  TFile *f1 = new TFile("../output/Data_ptGT30_TCHEM.root");
  TFile *f2 = new TFile("../output/MC_QCD-Pt15_ptGT30_TCHEM.root");

 char* xtitle = "TCHE discriminator";
 float binmin = 0., binmax = 25.;

 TH1F* h0= new TH1F("h0","",nbin,binmin,binmax);
 TH1F* h1= new TH1F("h1","",nbin,binmin,binmax);
 TH1F* h2= new TH1F("h2","",nbin,binmin,binmax);
 TH1F* h3= new TH1F("h3","",nbin,binmin,binmax);
 TH1F* g1= new TH1F("g1","",nbin,binmin,binmax);
 TH1F* g2= new TH1F("g2","",nbin,binmin,binmax);
 TH1F* g3= new TH1F("g3","",nbin,binmin,binmax);
      f1->cd();
 TH1F* h00 = (TH1F*)gROOT->FindObject("hData_Tagger");
      f2->cd();
 TH1F* g10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* g11 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_K0s");
 TH1F* g12 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Lam");
 TH1F* g20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* g30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");

       Float_t htot = 0, gtot = 0, g1tot = 0, g2tot = 0, g3tot = 0;
       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = h00->GetBinContent(i+nbin);
        h0->SetBinContent(i,y);
        htot += y;
        Float_t y  = g10->GetBinContent(i+nbin);
        Float_t y1 = g11->GetBinContent(i+nbin)*SFk0;
        Float_t y2 = g12->GetBinContent(i+nbin)*SFla;
        g1->SetBinContent(i,y+y1+y2);
        g1tot += y+y1+y2;
        Float_t y = g20->GetBinContent(i+nbin);
        g2->SetBinContent(i,y);
        g2tot += y;
        Float_t y = g30->GetBinContent(i+nbin);
        g3->SetBinContent(i,y);
        g3tot += y;
       }
       gtot = g1tot + g2tot + g3tot;
       Float_t LMC = g1tot / gtot;
       Float_t CMC = g2tot / gtot;
       Float_t BMC = g3tot / gtot;
       
       TObjArray *mc = new TObjArray(3); // MC histograms are put in this array
       mc->Add(g1);
       mc->Add(g2);
       mc->Add(g3);
       TFractionFitter* fitter = new TFractionFitter(h0, mc); // initialise
       fitter->Constrain(1,0.0,1.0);	 // constrain fraction 1 to be between 0 and 1
       fitter->Constrain(2,0.0,1.0);
       fitter->Constrain(3,0.0,1.0);
       fitter->SetRangeX(1,nbin);        // range for fit

       Int_t status = fitter->Fit();     // perform the fit
       cout << "###############" << endl;
       cout << " fit status: " << status << endl;
       cout << "###############" << endl;

       Double_t LFrac, CFrac, BFrac, LFracErr, CFracErr, BFracErr, Chi2, Ndf;
     if (status == 0) {                       // check on fit status
       TH1F* result = (TH1F*) fitter->GetPlot();
       fitter->GetResult(0, LFrac, LFracErr);
       fitter->GetResult(1, CFrac, CFracErr);
       fitter->GetResult(2, BFrac, BFracErr);
       Chi2 = fitter->GetChisquare();
       Ndf = fitter->GetNDF();

       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = g1->GetBinContent(i);
        g1->SetBinContent(i,y*LFrac*htot/g1tot);
        Float_t y = g2->GetBinContent(i);
        g2->SetBinContent(i,y*CFrac*htot/g2tot);
        Float_t y = g3->GetBinContent(i);
        g3->SetBinContent(i,y*BFrac*htot/g3tot);
       }
       h1->Add(g1,h1,1,0);
       h2->Add(g2,h1,1,1);
       h3->Add(g3,h2,1,1);

       h0->Draw("E");
       h0->SetLineColor(1);
       h0->SetLineStyle(1);
       h0->SetLineWidth(1);
       h0->SetMarkerStyle(20);
       h0->SetMarkerColor(1);
       h0->SetMarkerSize(0.5);
       h3->Draw("Hsame"); 
       h3->SetFillColor(kRed);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kGreen);
       h1->Draw("Hsame"); 
       h1->SetFillColor(kBlue);
       h0->SetTickLength(-0.02,"X");
       h0->SetLabelOffset(0.013,"X");
       h0->SetTickLength(-0.02,"Y");
       h0->SetLabelOffset(0.013,"Y");
       h0->GetXaxis()->SetLabelSize(0.07);
       h0->GetYaxis()->SetLabelSize(0.07);
       h0->GetXaxis()->SetTitleSize(0.07);
       h0->GetXaxis()->SetTitle(xtitle);
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h0->SetMaximum(3000000);  
       result->Draw("same");
       result->SetLineColor(1);
       result->SetLineStyle(1);
       result->SetLineWidth(2);
       h0->Draw("Esame"); 
     }

  cout << "##################################################" << endl;
  cout << " Chi2 / NDF = " << Chi2 << " / " << Ndf << endl;
  cout << "     b MC: " << BMC << " fit: " << BFrac << " +_ " << BFracErr << endl;
  cout << "     c MC: " << CMC << " fit: " << CFrac << " +_ " << CFracErr << endl;
  cout << " light MC: " << LMC << " fit: " << LFrac << " +_ " << LFracErr << endl;
  cout << "##################################################" << endl;

  Char_t * val = new Char_t[40]; // a pointer of type char

  TLegend* leg = new TLegend(0.33,0.73,0.93,0.93);
//     leg->SetHeader(" all jets");
    leg->AddEntry(h0," Data ","P");
  sprintf( val, " fit   #chi^{2} / ndf = %.2f / %.0f" , Chi2 , Ndf ) ;
    leg->AddEntry(result,val,"L");
  sprintf( val, "b-jet  %.3f #pm %.3f  (MC %.2f)" , BFrac , BFracErr , BMC) ;
    leg->AddEntry(h3,val,"F");
  sprintf( val, "c-jet  %.3f #pm %.3f  (MC %.2f)" , CFrac , CFracErr , CMC) ;
    leg->AddEntry(h2,val,"F");
  sprintf( val, "light  %.3f #pm %.3f  (MC %.2f)" , LFrac , LFracErr , LMC) ;
    leg->AddEntry(h1,val,"F");
    leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     pad2->cd();
  TFile *f1 = new TFile("../output/Data_ptGT30_TCHPM.root");
  TFile *f2 = new TFile("../output/MC_QCD-Pt15_ptGT30_TCHPM.root");

 char* xtitle = "TCHP discriminator";
 float binmin = 0., binmax = 25.;

 TH1F* h0= new TH1F("h0","",nbin,binmin,binmax);
 TH1F* h1= new TH1F("h1","",nbin,binmin,binmax);
 TH1F* h2= new TH1F("h2","",nbin,binmin,binmax);
 TH1F* h3= new TH1F("h3","",nbin,binmin,binmax);
 TH1F* g1= new TH1F("g1","",nbin,binmin,binmax);
 TH1F* g2= new TH1F("g2","",nbin,binmin,binmax);
 TH1F* g3= new TH1F("g3","",nbin,binmin,binmax);
      f1->cd();
 TH1F* h00 = (TH1F*)gROOT->FindObject("hData_Tagger");
      f2->cd();
 TH1F* g10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* g11 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_K0s");
 TH1F* g12 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Lam");
 TH1F* g20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* g30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");

       Float_t htot = 0, gtot = 0, g1tot = 0, g2tot = 0, g3tot = 0;
       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = h00->GetBinContent(i+nbin);
        h0->SetBinContent(i,y);
        htot += y;
        Float_t y  = g10->GetBinContent(i+nbin);
        Float_t y1 = g11->GetBinContent(i+nbin)*SFk0;
        Float_t y2 = g12->GetBinContent(i+nbin)*SFla;
        g1->SetBinContent(i,y+y1+y2);
        g1tot += y+y1+y2;
        Float_t y = g20->GetBinContent(i+nbin);
        g2->SetBinContent(i,y);
        g2tot += y;
        Float_t y = g30->GetBinContent(i+nbin);
        g3->SetBinContent(i,y);
        g3tot += y;
       }
       gtot = g1tot + g2tot + g3tot;
       Float_t LMC = g1tot / gtot;
       Float_t CMC = g2tot / gtot;
       Float_t BMC = g3tot / gtot;
       
       TObjArray *mc = new TObjArray(3); // MC histograms are put in this array
       mc->Add(g1);
       mc->Add(g2);
       mc->Add(g3);
       TFractionFitter* fitter = new TFractionFitter(h0, mc); // initialise
       fitter->Constrain(1,0.0,1.0);	 // constrain fraction 1 to be between 0 and 1
       fitter->Constrain(2,0.0,1.0);
       fitter->Constrain(3,0.0,1.0);
       fitter->SetRangeX(1,nbin);        // range for fit

       Int_t status = fitter->Fit();     // perform the fit
       cout << "###############" << endl;
       cout << " fit status: " << status << endl;
       cout << "###############" << endl;

       Double_t LFrac, CFrac, BFrac, LFracErr, CFracErr, BFracErr, Chi2, Ndf;
     if (status == 0) {                       // check on fit status
       TH1F* result = (TH1F*) fitter->GetPlot();
       fitter->GetResult(0, LFrac, LFracErr);
       fitter->GetResult(1, CFrac, CFracErr);
       fitter->GetResult(2, BFrac, BFracErr);
       Chi2 = fitter->GetChisquare();
       Ndf = fitter->GetNDF();

       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = g1->GetBinContent(i);
        g1->SetBinContent(i,y*LFrac*htot/g1tot);
        Float_t y = g2->GetBinContent(i);
        g2->SetBinContent(i,y*CFrac*htot/g2tot);
        Float_t y = g3->GetBinContent(i);
        g3->SetBinContent(i,y*BFrac*htot/g3tot);
       }
       h1->Add(g1,h1,1,0);
       h2->Add(g2,h1,1,1);
       h3->Add(g3,h2,1,1);

       h0->Draw("E");
       h0->SetLineColor(1);
       h0->SetLineStyle(1);
       h0->SetLineWidth(1);
       h0->SetMarkerStyle(20);
       h0->SetMarkerColor(1);
       h0->SetMarkerSize(0.5);
       h3->Draw("Hsame"); 
       h3->SetFillColor(kRed);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kGreen);
       h1->Draw("Hsame"); 
       h1->SetFillColor(kBlue);
       h0->SetTickLength(-0.02,"X");
       h0->SetLabelOffset(0.013,"X");
       h0->SetTickLength(-0.02,"Y");
       h0->SetLabelOffset(0.013,"Y");
       h0->GetXaxis()->SetLabelSize(0.07);
       h0->GetYaxis()->SetLabelSize(0.07);
       h0->GetXaxis()->SetTitleSize(0.07);
       h0->GetXaxis()->SetTitle(xtitle);
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h0->SetMaximum(3000000);  
       result->Draw("same");
       result->SetLineColor(1);
       result->SetLineStyle(1);
       result->SetLineWidth(2);
       h0->Draw("Esame"); 
     }

  cout << "##################################################" << endl;
  cout << " Chi2 / NDF = " << Chi2 << " / " << Ndf << endl;
  cout << "     b MC: " << BMC << " fit: " << BFrac << " +_ " << BFracErr << endl;
  cout << "     c MC: " << CMC << " fit: " << CFrac << " +_ " << CFracErr << endl;
  cout << " light MC: " << LMC << " fit: " << LFrac << " +_ " << LFracErr << endl;
  cout << "##################################################" << endl;

  Char_t * val = new Char_t[40]; // a pointer of type char

  TLegend* leg = new TLegend(0.33,0.73,0.93,0.93);
//     leg->SetHeader(" all jets");
    leg->AddEntry(h0," Data ","P");
  sprintf( val, " fit   #chi^{2} / ndf = %.2f / %.0f" , Chi2 , Ndf ) ;
    leg->AddEntry(result,val,"L");
  sprintf( val, "b-jet  %.3f #pm %.3f  (MC %.2f)" , BFrac , BFracErr , BMC) ;
    leg->AddEntry(h3,val,"F");
  sprintf( val, "c-jet  %.3f #pm %.3f  (MC %.2f)" , CFrac , CFracErr , CMC) ;
    leg->AddEntry(h2,val,"F");
  sprintf( val, "light  %.3f #pm %.3f  (MC %.2f)" , LFrac , LFracErr , LMC) ;
    leg->AddEntry(h1,val,"F");
    leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     pad3->cd();
  TFile *f1 = new TFile("../output/Data_ptGT30_SSVHEM.root");
  TFile *f2 = new TFile("../output/MC_QCD-Pt15_ptGT30_SSVHEM.root");

 char* xtitle = "SSVHE discriminator";
 float binmin = 0., binmax = 5.;
 int nbinew =40, binminew = 1.;

 TH1F* h0= new TH1F("h0","",nbinew,binminew,binmax);
 TH1F* h1= new TH1F("h1","",nbinew,binminew,binmax);
 TH1F* h2= new TH1F("h2","",nbinew,binminew,binmax);
 TH1F* h3= new TH1F("h3","",nbinew,binminew,binmax);
 TH1F* g1= new TH1F("g1","",nbinew,binminew,binmax);
 TH1F* g2= new TH1F("g2","",nbinew,binminew,binmax);
 TH1F* g3= new TH1F("g3","",nbinew,binminew,binmax);
      f1->cd();
 TH1F* h00 = (TH1F*)gROOT->FindObject("hData_Tagger");
      f2->cd();
 TH1F* g10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* g11 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_K0s");
 TH1F* g12 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Lam");
 TH1F* g20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* g30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");

       Float_t htot = 0, gtot = 0, g1tot = 0, g2tot = 0, g3tot = 0;
       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = h00->GetBinContent(i+nbin);
        h0->SetBinContent(i-10,y);
        htot += y;
        Float_t y  = g10->GetBinContent(i+nbin);
        Float_t y1 = g11->GetBinContent(i+nbin)*SFk0;
        Float_t y2 = g12->GetBinContent(i+nbin)*SFla;
        g1->SetBinContent(i-10,y+y1+y2);
        g1tot += y+y1+y2;
        Float_t y = g20->GetBinContent(i+nbin);
        g2->SetBinContent(i-10,y);
        g2tot += y;
        Float_t y = g30->GetBinContent(i+nbin);
        g3->SetBinContent(i-10,y);
        g3tot += y;
       }
       gtot = g1tot + g2tot + g3tot;
       Float_t LMC = g1tot / gtot;
       Float_t CMC = g2tot / gtot;
       Float_t BMC = g3tot / gtot;
       
       TObjArray *mc = new TObjArray(3); // MC histograms are put in this array
       mc->Add(g1);
       mc->Add(g2);
       mc->Add(g3);
       TFractionFitter* fitter = new TFractionFitter(h0, mc); // initialise
       fitter->Constrain(1,0.0,1.0);	 // constrain fraction 1 to be between 0 and 1
       fitter->Constrain(2,0.0,1.0);
       fitter->Constrain(3,0.0,1.0);
       fitter->SetRangeX(1,nbin);        // range for fit

       Int_t status = fitter->Fit();     // perform the fit
       cout << "###############" << endl;
       cout << " fit status: " << status << endl;
       cout << "###############" << endl;

       Double_t LFrac, CFrac, BFrac, LFracErr, CFracErr, BFracErr, Chi2, Ndf;
     if (status == 0) {                       // check on fit status
       TH1F* result = (TH1F*) fitter->GetPlot();
       fitter->GetResult(0, LFrac, LFracErr);
       fitter->GetResult(1, CFrac, CFracErr);
       fitter->GetResult(2, BFrac, BFracErr);
       Chi2 = fitter->GetChisquare();
       Ndf = fitter->GetNDF();

       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = g1->GetBinContent(i);
        g1->SetBinContent(i,y*LFrac*htot/g1tot);
        Float_t y = g2->GetBinContent(i);
        g2->SetBinContent(i,y*CFrac*htot/g2tot);
        Float_t y = g3->GetBinContent(i);
        g3->SetBinContent(i,y*BFrac*htot/g3tot);
       }
       h1->Add(g1,h1,1,0);
       h2->Add(g2,h1,1,1);
       h3->Add(g3,h2,1,1);

       h0->Draw("E");
       h0->SetLineColor(1);
       h0->SetLineStyle(1);
       h0->SetLineWidth(1);
       h0->SetMarkerStyle(20);
       h0->SetMarkerColor(1);
       h0->SetMarkerSize(0.5);
       h3->Draw("Hsame"); 
       h3->SetFillColor(kRed);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kGreen);
       h1->Draw("Hsame"); 
       h1->SetFillColor(kBlue);
       h0->SetTickLength(-0.02,"X");
       h0->SetLabelOffset(0.013,"X");
       h0->SetTickLength(-0.02,"Y");
       h0->SetLabelOffset(0.013,"Y");
       h0->GetXaxis()->SetLabelSize(0.07);
       h0->GetYaxis()->SetLabelSize(0.07);
       h0->GetXaxis()->SetTitleSize(0.07);
       h0->GetXaxis()->SetTitle(xtitle);
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h0->SetMaximum(3000000);  
       result->Draw("same");
       result->SetLineColor(1);
       result->SetLineStyle(1);
       result->SetLineWidth(2);
       h0->Draw("Esame"); 
     }

  cout << "##################################################" << endl;
  cout << " Chi2 / NDF = " << Chi2 << " / " << Ndf << endl;
  cout << "     b MC: " << BMC << " fit: " << BFrac << " +_ " << BFracErr << endl;
  cout << "     c MC: " << CMC << " fit: " << CFrac << " +_ " << CFracErr << endl;
  cout << " light MC: " << LMC << " fit: " << LFrac << " +_ " << LFracErr << endl;
  cout << "##################################################" << endl;

  Char_t * val = new Char_t[40]; // a pointer of type char

  TLegend* leg = new TLegend(0.33,0.73,0.93,0.93);
//     leg->SetHeader(" all jets");
    leg->AddEntry(h0," Data ","P");
  sprintf( val, " fit   #chi^{2} / ndf = %.2f / %.0f" , Chi2 , Ndf ) ;
    leg->AddEntry(result,val,"L");
  sprintf( val, "b-jet  %.3f #pm %.3f  (MC %.2f)" , BFrac , BFracErr , BMC) ;
    leg->AddEntry(h3,val,"F");
  sprintf( val, "c-jet  %.3f #pm %.3f  (MC %.2f)" , CFrac , CFracErr , CMC) ;
    leg->AddEntry(h2,val,"F");
  sprintf( val, "light  %.3f #pm %.3f  (MC %.2f)" , LFrac , LFracErr , LMC) ;
    leg->AddEntry(h1,val,"F");
    leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     pad4->cd();
  TFile *f1 = new TFile("../output/Data_ptGT30_SSVHPT.root");
  TFile *f2 = new TFile("../output/MC_QCD-Pt15_ptGT30_SSVHPT.root");

 char* xtitle = "SSVHP discriminator";
 float binmin = 0., binmax = 5.;
 int nbinew =40, binminew = 1.;

 TH1F* h0= new TH1F("h0","",nbinew,binminew,binmax);
 TH1F* h1= new TH1F("h1","",nbinew,binminew,binmax);
 TH1F* h2= new TH1F("h2","",nbinew,binminew,binmax);
 TH1F* h3= new TH1F("h3","",nbinew,binminew,binmax);
 TH1F* g1= new TH1F("g1","",nbinew,binminew,binmax);
 TH1F* g2= new TH1F("g2","",nbinew,binminew,binmax);
 TH1F* g3= new TH1F("g3","",nbinew,binminew,binmax);
      f1->cd();
 TH1F* h00 = (TH1F*)gROOT->FindObject("hData_Tagger");
      f2->cd();
 TH1F* g10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* g11 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_K0s");
 TH1F* g12 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger_Lam");
 TH1F* g20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* g30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");

       Float_t htot = 0, gtot = 0, g1tot = 0, g2tot = 0, g3tot = 0;
       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = h00->GetBinContent(i+nbin);
        h0->SetBinContent(i-10,y);
        htot += y;
        Float_t y  = g10->GetBinContent(i+nbin);
        Float_t y1 = g11->GetBinContent(i+nbin)*SFk0;
        Float_t y2 = g12->GetBinContent(i+nbin)*SFla;
        g1->SetBinContent(i-10,y+y1+y2);
        g1tot += y+y1+y2;
        Float_t y = g20->GetBinContent(i+nbin);
        g2->SetBinContent(i-10,y);
        g2tot += y;
        Float_t y = g30->GetBinContent(i+nbin);
        g3->SetBinContent(i-10,y);
        g3tot += y;
       }
       gtot = g1tot + g2tot + g3tot;
       Float_t LMC = g1tot / gtot;
       Float_t CMC = g2tot / gtot;
       Float_t BMC = g3tot / gtot;
       
       TObjArray *mc = new TObjArray(3); // MC histograms are put in this array
       mc->Add(g1);
       mc->Add(g2);
       mc->Add(g3);
       TFractionFitter* fitter = new TFractionFitter(h0, mc); // initialise
       fitter->Constrain(1,0.0,1.0);	 // constrain fraction 1 to be between 0 and 1
       fitter->Constrain(2,0.0,1.0);
       fitter->Constrain(3,0.0,1.0);
       fitter->SetRangeX(1,nbin);        // range for fit

       Int_t status = fitter->Fit();     // perform the fit
       cout << "###############" << endl;
       cout << " fit status: " << status << endl;
       cout << "###############" << endl;

       Double_t LFrac, CFrac, BFrac, LFracErr, CFracErr, BFracErr, Chi2, Ndf;
     if (status == 0) {                       // check on fit status
       TH1F* result = (TH1F*) fitter->GetPlot();
       fitter->GetResult(0, LFrac, LFracErr);
       fitter->GetResult(1, CFrac, CFracErr);
       fitter->GetResult(2, BFrac, BFracErr);
       Chi2 = fitter->GetChisquare();
       Ndf = fitter->GetNDF();

       for (Int_t i=1; i<=nbin; i++) {
        Float_t y = g1->GetBinContent(i);
        g1->SetBinContent(i,y*LFrac*htot/g1tot);
        Float_t y = g2->GetBinContent(i);
        g2->SetBinContent(i,y*CFrac*htot/g2tot);
        Float_t y = g3->GetBinContent(i);
        g3->SetBinContent(i,y*BFrac*htot/g3tot);
       }
       h1->Add(g1,h1,1,0);
       h2->Add(g2,h1,1,1);
       h3->Add(g3,h2,1,1);

       h0->Draw("E");
       h0->SetLineColor(1);
       h0->SetLineStyle(1);
       h0->SetLineWidth(1);
       h0->SetMarkerStyle(20);
       h0->SetMarkerColor(1);
       h0->SetMarkerSize(0.5);
       h3->Draw("Hsame"); 
       h3->SetFillColor(kRed);
       h2->Draw("Hsame"); 
       h2->SetFillColor(kGreen);
       h1->Draw("Hsame"); 
       h1->SetFillColor(kBlue);
       h0->SetTickLength(-0.02,"X");
       h0->SetLabelOffset(0.013,"X");
       h0->SetTickLength(-0.02,"Y");
       h0->SetLabelOffset(0.013,"Y");
       h0->GetXaxis()->SetLabelSize(0.07);
       h0->GetYaxis()->SetLabelSize(0.07);
       h0->GetXaxis()->SetTitleSize(0.07);
       h0->GetXaxis()->SetTitle(xtitle);
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h0->SetMaximum(3000000);  
       h0->SetMaximum(550);  
       result->Draw("same");
       result->SetLineColor(1);
       result->SetLineStyle(1);
       result->SetLineWidth(2);
       h0->Draw("Esame"); 
     }

  cout << "##################################################" << endl;
  cout << " Chi2 / NDF = " << Chi2 << " / " << Ndf << endl;
  cout << "     b MC: " << BMC << " fit: " << BFrac << " +_ " << BFracErr << endl;
  cout << "     c MC: " << CMC << " fit: " << CFrac << " +_ " << CFracErr << endl;
  cout << " light MC: " << LMC << " fit: " << LFrac << " +_ " << LFracErr << endl;
  cout << "##################################################" << endl;

  Char_t * val = new Char_t[40]; // a pointer of type char

  TLegend* leg = new TLegend(0.33,0.73,0.93,0.93);
//     leg->SetHeader(" all jets");
    leg->AddEntry(h0," Data ","P");
  sprintf( val, " fit   #chi^{2} / ndf = %.2f / %.0f" , Chi2 , Ndf ) ;
    leg->AddEntry(result,val,"L");
  sprintf( val, "b-jet  %.3f #pm %.3f  (MC %.2f)" , BFrac , BFracErr , BMC) ;
    leg->AddEntry(h3,val,"F");
  sprintf( val, "c-jet  %.3f #pm %.3f  (MC %.2f)" , CFrac , CFracErr , CMC) ;
    leg->AddEntry(h2,val,"F");
  sprintf( val, "light  %.3f #pm %.3f  (MC %.2f)" , LFrac , LFracErr , LMC) ;
    leg->AddEntry(h1,val,"F");
    leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
  c1->Update();
}

