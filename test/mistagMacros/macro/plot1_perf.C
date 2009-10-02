plot()
{
 TFile *f1  = new TFile("../output/QCDpt300_preprod_TCm.root");
 TFile *f11 = new TFile("../output/QCDpt300_preprod_TCt.root");
 TFile *f21 = new TFile("../output/QCDpt300_preprod_SVm.root");
 TFile *f31 = new TFile("../output/QCDpt300_preprod_JPm.root");
 TFile *f41 = new TFile("../output/QCDpt300_preprod_COm.root");

int stati=0;
bool fit= 0;
bool logy=1;

// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",200,10,700,730);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

// TPaveLabel *p01 = new TPaveLabel(0.1,0.94,0.78,0.98,
// //                       "MC qcd : p_{T}>20  #geq2 tracks","br");
//                       "QCD 80-120 Track Counting","br");
// p01->SetFillColor(0);
// p01->SetTextSize(0.8);
// p01->Draw();

pad1 = new TPad("pad1","This is pad1",0.04,0.04,0.92,0.92,21);

pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
if (logy) pad1->SetLogy();
pad1->SetGrid();

//$$ gStyle->SetOptDate(1);
gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleColor(29);
gStyle->SetTitleW(0.6);
gStyle->SetTitleH(0.05);
gStyle->SetOptStat(stati);

if (fit) {
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.1);
gStyle->SetOptFit(111);
} else {
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.2);
gStyle->SetOptFit(0);
}

// *****************************************************************************

  int nbin = 100;
  TVectorD cont(nbin);
  TVectorD beff(nbin);
  TVectorD ceff(nbin);
  TVectorD leff(nbin);
  float x1[100], ex1[100];
  float y1[100], ey1[100];
  float y2[100], ey2[100];

  pad1->cd();
//$$  TH2F* histo = new TH2F("histo","QCD 80-120 Track Counting"
  TH2F* histo = new TH2F("histo","",1000,0.,1.0,1000,0.0001,1.);
  histo->Draw(); 

/////////////////
// TC High Eff

  f1->cd();
  TH1F* h0 = (TH1F*)gROOT->FindObject("hBFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hBFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  beff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    beff(k) = cont(k) / h0->GetEntries();
    x1[k] = beff(k);
    ex1[k] = TMath::Sqrt(beff(k)*(1.-beff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hUDSFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  leff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    leff(k) = cont(k) / h0->GetEntries();
    y1[k] = leff(k);
    ey1[k] = TMath::Sqrt(leff(k)*(1.-leff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hCFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hCFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  ceff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    ceff(k) = cont(k) / h0->GetEntries();
    y2[k] = ceff(k);
    ey2[k] = TMath::Sqrt(ceff(k)*(1.-ceff(k))/h0->GetEntries());
  }  

  gr1 = new TGraphErrors(nbin,x1,y1,ex1,ey1);
  gr1->SetLineColor(4);
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(20);
  gr1->Draw("P"); 
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.04);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleColor(1);
  histo->GetYaxis()->SetTitleColor(1);
  histo->GetXaxis()->SetTitle("b efficiency");
//   histo->GetYaxis()->SetTitle("non-b efficiency");
  histo->GetYaxis()->SetTitle("udsg efficiency");
//   gr2 = new TGraphErrors(nbin,x1,y2,ex1,ey2);
//   histo->Draw("same"); 
//   gr2->SetMarkerColor(4);
//   gr2->SetMarkerStyle(24);
//   gr2->Draw("P"); 

/////////////////
// TC High Pur

  f11->cd();
  TH1F* h0 = (TH1F*)gROOT->FindObject("hBFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hBFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  beff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    beff(k) = cont(k) / h0->GetEntries();
    x1[k] = beff(k);
    ex1[k] = TMath::Sqrt(beff(k)*(1.-beff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hUDSFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  leff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    leff(k) = cont(k) / h0->GetEntries();
    y1[k] = leff(k);
    ey1[k] = TMath::Sqrt(leff(k)*(1.-leff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hCFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hCFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  ceff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    ceff(k) = cont(k) / h0->GetEntries();
    y2[k] = ceff(k);
    ey2[k] = TMath::Sqrt(ceff(k)*(1.-ceff(k))/h0->GetEntries());
  }  

  gr11 = new TGraphErrors(nbin,x1,y1,ex1,ey1);
  gr11->SetLineColor(7);
  gr11->SetMarkerColor(7);
  gr11->SetMarkerStyle(21);
  gr11->Draw("P"); 
//   gr12 = new TGraphErrors(nbin,x1,y2,ex1,ey2);
//   gr12->SetMarkerColor(7);
//   gr12->SetMarkerStyle(25);
//   gr12->Draw("P"); 

/////////////////
// Secondary Vertex

  f21->cd();
  TH1F* h0 = (TH1F*)gROOT->FindObject("hBFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hBFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  beff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    beff(k) = cont(k) / h0->GetEntries();
    x1[k] = beff(k);
    ex1[k] = TMath::Sqrt(beff(k)*(1.-beff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hUDSFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  leff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    leff(k) = cont(k) / h0->GetEntries();
    y1[k] = leff(k);
    ey1[k] = TMath::Sqrt(leff(k)*(1.-leff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hCFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hCFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  ceff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    ceff(k) = cont(k) / h0->GetEntries();
    y2[k] = ceff(k);
    ey2[k] = TMath::Sqrt(ceff(k)*(1.-ceff(k))/h0->GetEntries());
  }  

  gr21 = new TGraphErrors(nbin,x1,y1,ex1,ey1);
  gr21->SetLineColor(3);
  gr21->SetMarkerColor(3);
  gr21->SetMarkerStyle(29);
  gr21->Draw("P"); 
//   gr22 = new TGraphErrors(nbin,x1,y2,ex1,ey2);
//   gr22->SetMarkerColor(3);
//   gr22->SetMarkerStyle(30);
//   gr22->Draw("P"); 

/////////////////
// Jet Proba pos

  f31->cd();
  TH1F* h0 = (TH1F*)gROOT->FindObject("hBFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hBFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  beff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    beff(k) = cont(k) / h0->GetEntries();
    x1[k] = beff(k);
    ex1[k] = TMath::Sqrt(beff(k)*(1.-beff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hUDSFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  leff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    leff(k) = cont(k) / h0->GetEntries();
    if (y1[k] > 0.0001) y1[k] = leff(k);
    else y1[k] = 0.00001;
    ey1[k] = TMath::Sqrt(leff(k)*(1.-leff(k))/h0->GetEntries());
  }  

//   for (int k=1; k<nbin; k++) {
//    cout << "k " << k << "   beff " << beff(k) << "   leff " << leff(k) << endl;
//    cout << "k " << k << "   ex1   " << ex1[k]   << "   ey1   " << ey1[k]   << endl;
//   }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hCFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hCFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  ceff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    ceff(k) = cont(k) / h0->GetEntries();
    y2[k] = ceff(k);
    ey2[k] = TMath::Sqrt(ceff(k)*(1.-ceff(k))/h0->GetEntries());
  }  

  gr31 = new TGraphErrors(nbin,x1,y1,ex1,ey1);
  gr31->SetLineColor(2);
  gr31->SetMarkerColor(2);
  gr31->SetMarkerStyle(22);
  gr31->Draw("E"); 
//   gr32 = new TGraphErrors(nbin,x1,y2,ex1,ey2);
//   gr32->SetMarkerColor(2);
//   gr32->SetMarkerStyle(26);
//   gr32->Draw("P"); 

/////////////////
// Combined tag

  f41->cd();
  TH1F* h0 = (TH1F*)gROOT->FindObject("hBFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hBFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  beff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    beff(k) = cont(k) / h0->GetEntries();
    x1[k] = beff(k);
    ex1[k] = TMath::Sqrt(beff(k)*(1.-beff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hUDSFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hUDSFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  leff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    leff(k) = cont(k) / h0->GetEntries();
    y1[k] = leff(k);
    ey1[k] = TMath::Sqrt(leff(k)*(1.-leff(k))/h0->GetEntries());
  }  

  TH1F* h0 = (TH1F*)gROOT->FindObject("hCFlav_All_NTracks");
  TH1F* h1 = (TH1F*)gROOT->FindObject("hCFlav_PosTag");
  cont(0) = h1->GetEntries() - h1->GetBinContent(1) - h1->GetBinContent(0);
  ceff(0) = cont(0) / h0->GetEntries();
  for (int k=1; k<nbin; k++) {
    cont(k) = cont(k-1) - h1->GetBinContent(k+1);
    ceff(k) = cont(k) / h0->GetEntries();
    y2[k] = ceff(k);
    ey2[k] = TMath::Sqrt(ceff(k)*(1.-ceff(k))/h0->GetEntries());
  }  

  gr41 = new TGraphErrors(nbin,x1,y1,ex1,ey1);
  gr41->SetLineColor(1);
  gr41->SetMarkerColor(1);
  gr41->SetMarkerStyle(30);
  gr41->Draw("P"); 
//   gr42 = new TGraphErrors(nbin,x1,y2,ex1,ey2);
//   gr42->SetMarkerColor(6);
//   gr42->SetMarkerStyle(27);
//   gr42->Draw("P"); 

  gr31->Draw("Psame"); 

  TLegend* leg = new TLegend(0.12,0.70,0.35,0.89);
    leg->SetHeader(" QCDpt300_preprod");
    leg->AddEntry(gr1, "TCHE","P");
    leg->AddEntry(gr11,"TCHP","P");
    leg->AddEntry(gr31,"JP","P");
    leg->AddEntry(gr21,"SSV","P");
    leg->AddEntry(gr41,"CSV","P");
    leg->Draw();

    c1->Update();
}
