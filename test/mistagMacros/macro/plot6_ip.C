plot()
{
 TFile *f1 = new TFile("../output/QCDpt300_preprod_TCm.root");
 TFile *f2 = new TFile("../output/QCDpt300_preprod_TCt.root");
 TFile *f3 = new TFile("../output/QCDpt300_preprod_JPm.root");
 TFile *f4 = new TFile("../output/QCDpt300_preprod_SVm.root");
 TFile *f5 = new TFile("../output/QCDpt300_preprod_COm.root");
 TFile *f6 = new TFile("../output/QCDpt300_preprod_MUm.root");

// *****************************************************************************

 Int_t stati=0;
 Bool_t  fit=0;
 Bool_t logy=1;

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
gStyle->SetTitleW(0.3);
gStyle->SetTitleH(0.08);
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
 TH1F* h0= new TH1F("h0","TCHE",100,-50.,50.);
 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* g1= new TH1F("g1","",100,-50.,50.);
 TH1F* g2= new TH1F("g2","",100,-50.,50.);
 TH1F* g3= new TH1F("g3","",100,-50.,50.);
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Tagger");
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");
       h0->Add(h00,h0,1,0);
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h2,1,1);
       for (Int_t i=1; i<51; i++) {
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
       h0->GetXaxis()->SetLabelSize(0.08);
       h0->GetYaxis()->SetLabelSize(0.08);
       h0->GetXaxis()->SetTitleSize(0.08);
       h0->GetXaxis()->SetTitle("IP/#sigma");
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h1->SetMaximum(3000000);  

  TLegend* leg = new TLegend(0.65,0.65,0.90,0.93);
//     leg->SetHeader(" all jets");
    leg->AddEntry(h3,"  b","F");
    leg->AddEntry(h2,"  c","F");
    leg->AddEntry(h1," udsg","F");
    leg->Draw();

    pad2->cd();
      f2->cd();
 TH1F* h0= new TH1F("h0","TCHP",100,-50.,50.);
 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* g1= new TH1F("g1","",100,-50.,50.);
 TH1F* g2= new TH1F("g2","",100,-50.,50.);
 TH1F* g3= new TH1F("g3","",100,-50.,50.);
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Tagger");
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");
       h0->Add(h00,h0,1,0);
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h2,1,1);
       for (Int_t i=1; i<51; i++) {
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
       h0->GetXaxis()->SetLabelSize(0.08);
       h0->GetYaxis()->SetLabelSize(0.08);
       h0->GetXaxis()->SetTitleSize(0.08);
       h0->GetXaxis()->SetTitle("IP/#sigma");
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h1->SetMaximum(3000000);  

    pad3->cd();
      f3->cd();
 TH1F* h0= new TH1F("h0","JP",100,-50.,50.);
 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* g1= new TH1F("g1","",100,-50.,50.);
 TH1F* g2= new TH1F("g2","",100,-50.,50.);
 TH1F* g3= new TH1F("g3","",100,-50.,50.);
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Tagger");
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");
       h0->Add(h00,h0,1,0);
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h2,1,1);
       for (Int_t i=1; i<51; i++) {
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
       h0->GetXaxis()->SetLabelSize(0.08);
       h0->GetYaxis()->SetLabelSize(0.08);
       h0->GetXaxis()->SetTitleSize(0.08);
       h0->GetXaxis()->SetTitle("Discri #times 20");
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h1->SetMaximum(3000000);  

    pad4->cd();
      f4->cd();
 TH1F* h0= new TH1F("h0","SSV",100,-50.,50.);
 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* g1= new TH1F("g1","",100,-50.,50.);
 TH1F* g2= new TH1F("g2","",100,-50.,50.);
 TH1F* g3= new TH1F("g3","",100,-50.,50.);
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Tagger");
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");
       h0->Add(h00,h0,1,0);
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h2,1,1);
       for (Int_t i=1; i<51; i++) {
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
       h0->GetXaxis()->SetLabelSize(0.08);
       h0->GetYaxis()->SetLabelSize(0.08);
       h0->GetXaxis()->SetTitleSize(0.08);
       h0->GetXaxis()->SetTitle("(Discri - 1) #times 10");
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h1->SetMaximum(3000000);  

    pad5->cd();
      f5->cd();
 TH1F* h0= new TH1F("h0","CSV",100,-50.,50.);
 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* g1= new TH1F("g1","",100,-50.,50.);
 TH1F* g2= new TH1F("g2","",100,-50.,50.);
 TH1F* g3= new TH1F("g3","",100,-50.,50.);
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Tagger");
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");
       h0->Add(h00,h0,1,0);
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h2,1,1);
       for (Int_t i=1; i<51; i++) {
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
       h0->GetXaxis()->SetLabelSize(0.08);
       h0->GetYaxis()->SetLabelSize(0.08);
       h0->GetXaxis()->SetTitleSize(0.08);
       h0->GetXaxis()->SetTitle("Discri #times 50");
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h1->SetMaximum(3000000);  

    pad6->cd();
      f6->cd();
 TH1F* h0= new TH1F("h0","Muon NN",100,-50.,50.);
 TH1F* h1= new TH1F("h1","",100,-50.,50.);
 TH1F* h2= new TH1F("h2","",100,-50.,50.);
 TH1F* h3= new TH1F("h3","",100,-50.,50.);
 TH1F* g1= new TH1F("g1","",100,-50.,50.);
 TH1F* g2= new TH1F("g2","",100,-50.,50.);
 TH1F* g3= new TH1F("g3","",100,-50.,50.);
 TH1F* h00 = (TH1F*)gROOT->FindObject("hAllFlav_Tagger");
 TH1F* h10 = (TH1F*)gROOT->FindObject("hLightFlav_Tagger");
 TH1F* h20 = (TH1F*)gROOT->FindObject("hCFlav_Tagger");
 TH1F* h30 = (TH1F*)gROOT->FindObject("hBFlav_Tagger");
       h0->Add(h00,h0,1,0);
       h1->Add(h10,h1,1,0);
       h2->Add(h20,h1,1,1);
       h3->Add(h30,h2,1,1);
       for (Int_t i=1; i<51; i++) {
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
       h0->GetXaxis()->SetLabelSize(0.08);
       h0->GetYaxis()->SetLabelSize(0.08);
       h0->GetXaxis()->SetTitleSize(0.08);
       h0->GetXaxis()->SetTitle("Discri #times 50");
       h0->GetXaxis()->SetTitleColor(1);
       h0->GetXaxis()->SetNdivisions(509);
       h0->GetYaxis()->SetNdivisions(509);
//        h0->SetMinimum(20); h1->SetMaximum(3000000);  

  c1->Update();
}

