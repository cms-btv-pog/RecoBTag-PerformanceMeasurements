plot()
{
// Primary MC sample
 TFile *f1 = new TFile("../output/QCDpt80_JPm.root");
// MC with Track History
 TFile *f2 = new TFile("../output/QCDpt300_preprod_JPm.root");

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

pad1 = new TPad("pad1","This is pad1",0.04,0.66,0.96,0.94,21);
pad2 = new TPad("pad2","This is pad2",0.04,0.35,0.96,0.63,21);
pad3 = new TPad("pad3","This is pad3",0.04,0.04,0.96,0.32,21);

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

pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(0);
   pad3->SetTopMargin(0.05);
   pad3->SetBottomMargin(0.15);
   pad3->SetRightMargin(0.05);
   pad3->SetLeftMargin(0.15);

//$$ gStyle->SetOptDate(1);
gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleColor(29);
gStyle->SetTitleW(0.4);
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

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Initialise Histograms

 TF1 *fun = new TF1("fun","[0]+[1]*x+[2]*x*x",30.,230.); // Loose Medium
//$$ TF1 *fun = new TF1("fun","([0]+[1]*x)/(1.+[2]*x*x)",30.,230.); // Tight
 TF1 *gun = new TF1("gun","[0]",80.,120.);
 TF1 *hun = new TF1("hun","[0]",40.,230.);
//  TF1 *hun = new TF1("hun","[0]+[1]*x*x+[2]*x*x*x",30.,230.);
      hun->SetLineColor(1);

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

//$$
 TH1F* h4= new TH1F("h4","",20,30.,230.);  // light pos / light
//  TH1F* h4= new TH1F("h4","Track Count Tight : no veto",20,30.,230.);  // light pos / light
//  TH1F* h4= new TH1F("h4","Jet Proba pos Medium : no veto",20,30.,230.);  // light pos / light
//  TH1F* h4= new TH1F("h4","Jet Proba all Medium : no veto",20,30.,230.);  // light pos / light
//  TH1F* h4= new TH1F("h4","Second Vertex Medium : no veto",20,30.,230.);  // light pos / light
       h4->Divide(h29,g09,1,1,"B"); 
//$$
 TH1F* h5= new TH1F("h5","",20,30.,230.);  // all neg with veto / all with veto
       h5->Divide(h10,h00,1,1,"B");

//  TH1F* h6= new TH1F("h6","Rl  =  mistag / neg.tag",20,30.,230.);  // (light pos/light) / (all neg/all)
 TH1F* h6= new TH1F("h6","",20,30.,230.);  // (light pos/light) / (all neg/all)
       h6->Divide(h4,h5,1,1);

 TH1F* h8= new TH1F("h8","",20,30.,230.);  // uds pos / uds
       h8->Divide(h21,g01,1,1,"B");

 TH1F* h9= new TH1F("h9","",20,30.,230.);  // g pos / g
       h9->Divide(h24,g04,1,1,"B");

     pad1->cd();
       h4->Draw("E"); 
       h4->SetLineColor(1);
       h4->SetMarkerStyle(20);
       h4->SetMarkerColor(1);
       h4->SetMarkerSize(1.);
       h5->Draw("Esame");
       h5->SetLineColor(6);
       h5->SetMarkerStyle(21);
       h5->SetMarkerColor(6);
       h5->SetMarkerSize(1.1);
//        h8->Draw("Esame");
//        h8->SetLineColor(3);
//        h8->SetMarkerStyle(23);
//        h8->SetMarkerColor(3);
//        h8->SetMarkerSize(1.);
//        h9->Draw("Esame");
//        h9->SetLineColor(3);
//        h9->SetMarkerStyle(22);
//        h9->SetMarkerColor(3);
//        h9->SetMarkerSize(1.);
//        h4->Draw("Esame"); 
//        h5->Draw("Esame");
//$$       
//$$       h5->Fit("gun","rvee"); 
//$$       
       h4->GetXaxis()->SetLabelSize(0.08);
       h4->GetYaxis()->SetLabelSize(0.08);
       h4->GetXaxis()->SetTitleSize(0.07);
       h4->GetXaxis()->SetTitle("p_{T}(jet) (GeV)");
       h4->GetXaxis()->SetTitleColor(1);
//        h4->SetMinimum(0.01); h4->SetMaximum(1.00); // Loose
       h4->SetMinimum(0.001); h4->SetMaximum(0.2); // Medium
//        h4->SetMinimum(0.0001); h4->SetMaximum(0.02); // Tight
       h4->GetXaxis()->SetNdivisions(509);
       h4->GetYaxis()->SetNdivisions(509);
       h4->Draw("Esame"); 

  TLegend* leg = new TLegend(0.20,0.73,0.40,0.93);
//   TLegend* leg = new TLegend(0.79,0.80,0.99,1.0);
//     leg->SetHeader("taggable");
//     leg->AddEntry(h9,"  g  mistag","P");
//     leg->AddEntry(h8,"uds mistag","P");
    leg->AddEntry(h4,"light mistag","P");
    leg->AddEntry(h5,"negative tag","P");
    leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Flavour fraction

 TH1F* k00= new TH1F("k00","",20,30.,230.);  // g/uds ratio
 TH1F* k01= new TH1F("k01","",20,30.,230.);  // uds-fraction
 TH1F* k02= new TH1F("k02","",20,30.,230.);  // c-fraction
 TH1F* k03= new TH1F("k03","",20,30.,230.);  // b-fraction
 TH1F* k04= new TH1F("k04","flavour fraction",20,30.,230.);  // g-fraction
 TH1F* k11= new TH1F("k11","",20,30.,230.);  // uds-fraction neg
 TH1F* k12= new TH1F("k12","flavour fraction",20,30.,230.);  // c-fraction neg
 TH1F* k13= new TH1F("k13","",20,30.,230.);  // b-fraction neg
 TH1F* k14= new TH1F("k14","",20,30.,230.);  // g-fraction neg
 TH1F* k21= new TH1F("k21","",20,30.,230.);  // uds-fraction pos
 TH1F* k24= new TH1F("k24","",20,30.,230.);  // g-fraction pos

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

//        k04->Draw("E"); 
//        k04->SetLineColor(3);
//        k04->SetMarkerStyle(26);
//        k04->SetMarkerColor(3);
//        k04->SetMarkerSize(1.);
//        k14->Draw("Esame"); 
//        k14->SetLineColor(3);
//        k14->SetMarkerStyle(22);
//        k14->SetMarkerColor(3);
//        k14->SetMarkerSize(1.);
//        k04->GetXaxis()->SetLabelSize(0.08);
//        k04->GetYaxis()->SetLabelSize(0.08);
//        k04->GetXaxis()->SetTitleSize(0.07);
//        k04->GetXaxis()->SetTitle("p_{T}(jet) (GeV)");
//        k04->GetXaxis()->SetTitleColor(1);
//        k04->GetXaxis()->SetNdivisions(509);
//        k04->GetYaxis()->SetNdivisions(509);
//        k04->SetMinimum(0.005); k04->SetMaximum(2.0);  
// //$$
// //$$       k04->Fit("hun","rvee"); 
// //$$
//        
//        k03->Draw("Esame"); 
//        k03->SetLineColor(2);
//        k03->SetMarkerStyle(24);
//        k03->SetMarkerColor(2);
//        k03->SetMarkerSize(1.);
// //$$
// //$$       k03->Fit("hun","rvee"); 
// //$$
//        k13->Draw("Esame"); 
//        k13->SetLineColor(2);
//        k13->SetMarkerStyle(20);
//        k13->SetMarkerColor(2);
//        k13->SetMarkerSize(1.);
// 
//        k02->Draw("Esame"); 
//        k02->SetLineColor(4);
//        k02->SetMarkerStyle(25);
//        k02->SetMarkerColor(4);
//        k02->SetMarkerSize(1.);
// //$$
// //$$       k02->Fit("hun","rvee"); 
// //$$
//        k12->Draw("Esame"); 
//        k12->SetLineColor(4);
//        k12->SetMarkerStyle(21);
//        k12->SetMarkerColor(4);
//        k12->SetMarkerSize(1.);
// 
//   TLegend* leg = new TLegend(0.82,0.80,0.99,1.0);
// //     leg->SetHeader("taggable");
//     leg->AddEntry(k14,"    g","P");
//     leg->AddEntry(k12,"    c","P");
//     leg->AddEntry(k13,"    b","P");
//     leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Systematics

// systematics due to flavour fraction

 TH1F* bb= new TH1F("bb","",20,30.,230.);
       for (int i=0; i<30; i++) {
        float yneg = k13->GetBinContent(i);
	float y = k13->GetBinContent(i) - k03->GetBinContent(i);
	y = Xbb * TMath::Abs( y );
	float ey = Xbb * k13->GetBinError(i);
        bb->SetBinContent(i,y);
        bb->SetBinError(i,ey);
       }
 TH1F* cc= new TH1F("cc","",20,30.,230.);
       for (int i=0; i<30; i++) {
        float yneg = k12->GetBinContent(i);
	float y = k12->GetBinContent(i) - k02->GetBinContent(i);
	y = Xcc * TMath::Abs( y );
	float ey = Xcc * k12->GetBinError(i);
        cc->SetBinContent(i,y);
        cc->SetBinError(i,ey);
       }

 TH1F* gg= new TH1F("gg","",20,30.,230.);
       for (int i=0; i<30; i++) {
        float frac = k00->GetBinContent(i);
        float lneg = k11->GetBinContent(i);
        float lpos = k21->GetBinContent(i);
        float gneg = k14->GetBinContent(i);
        float gpos = k24->GetBinContent(i);
	float yl = frac * (lpos - lneg);
	float yg = gpos - gneg;
	y = Xgg * TMath::Abs( yg - yl );
// cout << i << " frac " << frac << " lneg " << lneg << " lpos " << lpos 
//                               << " gneg " << gneg << " gpos " << gpos
// 			      << " yneg " << yneg << " ypos " << ypos 
// 			      << " y " << y << endl;
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

 TH1F* h1= new TH1F("h1","",20,30.,230.); // all pos / all neg
       h1->Divide(h20,g10,1,1);
 TH1F* h2= new TH1F("h2","",20,30.,230.); // light neg / light pos
       h2->Divide(g19,h29,1,1);
 TH1F* sf= new TH1F("sf","",20,30.,230.); // sign flip
       sf->Add(h1,h2,Xip,-Xip);

// systematics due to V0

 TH1F* h71= new TH1F("h71","",20,30.,230.); // K0 neg / all neg
       h71->Divide(h110,i10,1,1,"B");
 TH1F* h72= new TH1F("h72","",20,30.,230.); // light K0 pos / light pos
       h72->Divide(h129,i29,1,1,"B");
 TH1F* h73= new TH1F("h73","",20,30.,230.); // light K0 / light
       h73->Divide(h109,i09,1,1,"B");
 TH1F* h74= new TH1F("h74","",20,30.,230.); // K0 / all
       h74->Divide(h100,i00,1,1,"B");

 TH1F* h81= new TH1F("h81","",20,30.,230.); // Bad neg / all neg
       h81->Divide(h210,i10,SFbad,1,"B");
 TH1F* h82= new TH1F("h82","",20,30.,230.); // light Bad pos / light pos
       h82->Divide(h229,i29,SFbad,1,"B");
 TH1F* h83= new TH1F("h83","",20,30.,230.); // light Bad / light
       h83->Divide(h209,i09,SFbad,1,"B");
 TH1F* h84= new TH1F("h84","",20,30.,230.); // Bad / all
       h84->Divide(h200,i00,SFbad,1,"B");

 TH1F* h91= new TH1F("h91","",20,30.,230.); // gamma neg / all neg
       h91->Divide(h310,i10,1,1,"B");
 TH1F* h92= new TH1F("h92","",20,30.,230.); // light gamma pos / light pos
       h92->Divide(h329,i29,1,1,"B");
 TH1F* h93= new TH1F("h93","",20,30.,230.); // light gamma / light
       h93->Divide(h309,i09,1,1,"B");
 TH1F* h94= new TH1F("h94","",20,30.,230.); // gamma / all
       h94->Divide(h300,i00,1,1,"B");

 TH1F* ka = new TH1F("ka","",20,30.,230.);
 TH1F* bad = new TH1F("bad","",20,30.,230.);
 TH1F* ga = new TH1F("ga","",20,30.,230.);
 TH1F* v0 = new TH1F("v0","",20,30.,230.);

       for (int i=0; i<30; i++) {
	float y    = h72->GetBinContent(i) - h71->GetBinContent(i);
	             - h73->GetBinContent(i) + h74->GetBinContent(i);
	y = Xv * TMath::Abs( y );
	float ey = h71->GetBinError(i)*h71->GetBinError(i) 
	           + h72->GetBinError(i)*h72->GetBinError(i);
	ey = Xv * TMath::Sqrt( ey );
        ka->SetBinContent(i,y);
        ka->SetBinError(i,ey);
       }
       for (int i=0; i<30; i++) {
	float y    = h82->GetBinContent(i) - h81->GetBinContent(i);
	             - h83->GetBinContent(i) + h84->GetBinContent(i);
	y = Xbad * TMath::Abs( y );
	float ey = h81->GetBinError(i)*h81->GetBinError(i) 
	           + h82->GetBinError(i)*h82->GetBinError(i);
	ey = Xbad * TMath::Sqrt( ey );
        bad->SetBinContent(i,y);
        bad->SetBinError(i,ey);
       }
       for (int i=0; i<30; i++) {
	float y    = h92->GetBinContent(i) - h91->GetBinContent(i);
	             - h93->GetBinContent(i) + h94->GetBinContent(i);
	y = Xv * TMath::Abs( y );
	float ey = h91->GetBinError(i)*h91->GetBinError(i) 
	           + h92->GetBinError(i)*h92->GetBinError(i);
	ey = Xv * TMath::Sqrt( ey );
        ga->SetBinContent(i,y);
        ga->SetBinError(i,ey);
       }
       for (int i=0; i<30; i++) {
	float y = ka->GetBinContent(i)*ka->GetBinContent(i) 
	          + ga->GetBinContent(i)*ga->GetBinContent(i);
	y = TMath::Sqrt( y );
	float ey = ka->GetBinError(i)*ka->GetBinError(i) 
	           + ga->GetBinError(i)*ga->GetBinError(i);
	ey = TMath::Sqrt( ey );
        v0->SetBinContent(i,y);
        v0->SetBinError(i,ey);
       }

// cout << "############################################################" << endl; 
// cout << "############################################################" << endl; 
//  TH1F* tot= new TH1F("tot","Systematics",20,30.,230.);
 TH1F* tot= new TH1F("tot","",20,30.,230.);
       for (int i=0; i<30; i++) {
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
// if (i == 7)  cout << " 85 GeV: b = " << ybb << "  c = " << ycc 
//                          << "  g = " << ygg << "  stat = " << ey << endl; 
// if (i == 8)  cout << " 95 GeV: b = " << ybb << "  c = " << ycc 
//                          << "  g = " << ygg << "  stat = " << ey << endl; 
// if (i == 9)  cout << "105 GeV: b = " << ybb << "  c = " << ycc 
//                          << "  g = " << ygg << "  stat = " << ey << endl; 
// if (i == 10) cout << "115 GeV: b = " << ybb << "  c = " << ycc 
//                          << "  g = " << ygg << "  stat = " << ey << endl; 
       }
// cout << "############################################################" << endl; 
// cout << "############################################################" << endl; 

     pad2->cd();

       tot->Draw("E"); 
       tot->SetLineColor(1);
       tot->SetMarkerStyle(23);
       tot->SetMarkerColor(1);
       tot->SetMarkerSize(1.);
       tot->GetXaxis()->SetLabelSize(0.08);
       tot->GetYaxis()->SetLabelSize(0.08);
       tot->GetXaxis()->SetTitleSize(0.07);
       tot->GetXaxis()->SetTitle("p_{T}(jet) (GeV)");
       tot->GetXaxis()->SetTitleColor(1);
       tot->GetXaxis()->SetNdivisions(509);
       tot->GetYaxis()->SetNdivisions(509);
//	  fun->SetLineColor(1);
//  	  tot->Fit("fun",""); 
//	  tot->Fit("hun","rvee"); 
//        tot->SetMinimum(0.); tot->SetMaximum(0.20); // Loose
       tot->SetMinimum(0.); tot->SetMaximum(0.60); // Medium
//        tot->SetMinimum(0.); tot->SetMaximum(0.80); // Tight
       
       bb->Draw("Esame"); 
       bb->SetLineColor(2);
       bb->SetMarkerStyle(20);
       bb->SetMarkerColor(2);
       bb->SetMarkerSize(1.);
//$$       bb->Fit("hun","rvee"); 

       cc->Draw("Esame"); 
       cc->SetLineColor(4);
       cc->SetMarkerStyle(21);
       cc->SetMarkerColor(4);
       cc->SetMarkerSize(1.);
//$$       cc->Fit("hun","rvee"); 

       gg->Draw("Esame"); 
       gg->SetLineColor(3);
       gg->SetMarkerStyle(22);
       gg->SetMarkerColor(3);
       gg->SetMarkerSize(1.);
//$$       gg->Fit("hun","rvee"); 

       v0->Draw("Esame"); 
       v0->SetLineColor(1);
       v0->SetMarkerStyle(24);
       v0->SetMarkerColor(1);
       v0->SetMarkerSize(1.);

       bad->Draw("Esame"); 
       bad->SetLineColor(1);
       bad->SetMarkerStyle(25);
       bad->SetMarkerColor(1);
       bad->SetMarkerSize(1.);

       sf->Draw("Esame"); 
       sf->SetLineColor(1);
       sf->SetMarkerStyle(30);
       sf->SetMarkerColor(1);
       sf->SetMarkerSize(1.);
 
       tot->Draw("Esame"); 

//   TLegend* leg = new TLegend(0.74,0.43,0.94,0.93);
//     leg->AddEntry(bb," b #pm 50% ","P");
//     leg->AddEntry(cc," c #pm 50% ","P");
//     leg->AddEntry(gg," g #pm 20% ","P");
//     leg->AddEntry(v0,"V^{0}+Int #pm 20%","P");
//     leg->AddEntry(bad," bad#times2 #pm 30%","P");
//     leg->AddEntry(sf," sign flip ","P");
//     leg->Draw();

  TLegend* leg = new TLegend(0.20,0.83,0.34,0.93);
   leg->AddEntry(tot,"all syst.","P");
  leg->Draw();

  TLegend* leg = new TLegend(0.36,0.68,0.50,0.93);
    leg->AddEntry(bb," b #pm 50% ","P");
    leg->AddEntry(cc," c #pm 50% ","P");
    leg->AddEntry(gg," g #pm 20% ","P");
    leg->Draw();

  TLegend* leg = new TLegend(0.52,0.68,0.72,0.93);
    leg->AddEntry(v0,"V^{0}+Int #pm 20%","P");
    leg->AddEntry(bad," bad+fake #pm 30%","P");
    leg->AddEntry(sf," sign flip ","P");
    leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Rlight

 TF1 *fun = new TF1("fun","[0]+[1]*x+[2]*x*x",30.,230.);
      fun->SetLineColor(1);
      fun->SetLineStyle(1);
 TF1 *gun = new TF1("gun","[0]+[1]*x+[2]*x*x",30.,230.);
      gun->SetLineColor(1);
      gun->SetLineStyle(2);
 TF1 *hun = new TF1("hun","[0]+[1]*x+[2]*x*x",30.,230.);
      hun->SetLineColor(1);
      hun->SetLineStyle(2);

 TH1F* u6= new TH1F("u6","",20,30.,230.);
 TH1F* d6= new TH1F("d6","",20,30.,230.);
       for (int i=0; i<30; i++) {
	float y = h6->GetBinContent(i);
	float ey = h6->GetBinError(i);
	float syst = tot->GetBinContent(i);
	float yu = y + TMath::Sqrt( ey*ey + syst*syst*y*y );
	float yd = y - TMath::Sqrt( ey*ey + syst*syst*y*y );
        u6->SetBinContent(i,yu);
        u6->SetBinError(i,ey);
        d6->SetBinContent(i,yd);
        d6->SetBinError(i,ey);
       }

     pad3->cd();
       h6->Draw("E"); 
       h6->SetLineColor(1);
       h6->SetMarkerStyle(20);
       h6->SetMarkerColor(1);
       h6->SetMarkerSize(1.);
       h6->GetXaxis()->SetLabelSize(0.08);
       h6->GetYaxis()->SetLabelSize(0.08);
       h6->GetXaxis()->SetTitleSize(0.07);
       h6->GetXaxis()->SetTitle("p_{T}(jet) (GeV)");
       h6->GetXaxis()->SetTitleColor(1);
       h6->GetXaxis()->SetNdivisions(509);
       h6->GetYaxis()->SetNdivisions(509);
       h6->Fit("fun","rvee"); 
       u6->Fit("gun","rvee0");
       gun->Draw("same"); 
       d6->Fit("hun","rvee0");
       hun->Draw("same"); 
//        h6->SetMinimum(0.5); h6->SetMaximum(2.5); // Loose
//        h6->SetMinimum(0.); h6->SetMaximum(5.); // Medium or Tight
       h6->SetMinimum(0.); h6->SetMaximum(2.); // Medium
//        h6->SetMinimum(0.); h6->SetMaximum(2.5); // Tight
//        h6->SetMinimum(0.); h6->SetMaximum(10.); // Tight veto

  TLegend* leg = new TLegend(0.20,0.73,0.45,0.93);
//   leg->SetHeader("Rl = mistag/neg.tag");
   leg->AddEntry(h6,"mistag/neg.tag","P");
   leg->AddEntry(gun,"#pm systematics","L");
  leg->Draw();

//  TH1F* h89= new TH1F("h89","",20,30.,230.);  // uds tag / g tag
//        h89->Divide(h8,h9,1,1);
//        h89->Draw("Esame");
//        h89->SetLineColor(3);
//        h89->SetMarkerStyle(23);
//        h89->SetMarkerColor(3);
//        h89->SetMarkerSize(1.);
// //$$
// //$$       h89->Fit("hun","rvee"); 
// //$$

//   TLegend* leg = new TLegend(0.42,0.75,0.92,0.90);
// //     leg->SetHeader("taggable");
//     leg->AddEntry(h6,"R_{light}","P");
//     leg->AddEntry(h89,"uds / g mistags","P");
//     leg->Draw();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Final Histograms

  TH1F* hMistag_JetPt     = new TH1F("hMistag_JetPt","pt(jet)",20,30.,230.);
  TH1F* hNegTag_JetPt     = new TH1F("hNegTag_JetPt","pt(jet)",20,30.,230.);
  TH1F* hRlight_JetPt     = new TH1F("hRlight_JetPt","pt(jet)",20,30.,230.);
  TH1F* hRlightSyst_JetPt = new TH1F("hRlightSyst_JetPt","pt(jet)",20,30.,230.);

  hMistag_JetPt -> Add(h4,h4,1,0);
  hNegTag_JetPt -> Add(h5,h5,1,0);
  hRlight_JetPt -> Add(h6,h6,1,0);
  hRlightSyst_JetPt -> Add(tot,tot,1,0);

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  c1->Update();
}

