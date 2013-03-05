{
   gROOT->ProcessLine(".L CommPlotProducer4ttbar.C++");
   
   // Lumi = 12.187/fb
   float lumi = 12187; // /pb

   // initialisation
   TFile* file_data = new TFile("/opt/sbg/data/data3/cms/ccollard/files/btag_ntuples/prod_ttbar/mueg_run2012_abc.root");
   TChain* tree_data = (TChain*) file_data->Get("btagana/ttree");
   CommPlotProducer4ttbar m_data(tree_data);

   double n_ttbar=0;
   double n_dy1=0;
   double n_dy2=0;
   double n_st1=0;
   double n_st2=0;

   TFile* file_ttbar = new TFile("/opt/sbg/data/data3/cms/ccollard/files/btag_ntuples/prod_ttbar/ttbar.root");
   TH1F* inputNumber = (TH1F*)file_ttbar->Get("ttbarselectionproducer/hcheck_cutflow");
   n_ttbar = inputNumber->GetBinContent(1);
   TChain* tree_ttbar = (TChain*) file_ttbar->Get("btagana/ttree");
   CommPlotProducer4ttbar m_ttbar(tree_ttbar);

   TFile* file_dy1 = new TFile("/opt/sbg/data/data3/cms/ccollard/files/btag_ntuples/prod_ttbar/dy_m10to50.root");
   TH1F* inputNumber2 = (TH1F*)file_dy1->Get("ttbarselectionproducer/hcheck_cutflow");
   n_dy1 = inputNumber2->GetBinContent(1);
   TChain* tree_dy1 = (TChain*) file_dy1->Get("btagana/ttree");
   CommPlotProducer4ttbar m_dy1(tree_dy1);

   TFile* file_dy2 = new TFile("/opt/sbg/data/data3/cms/ccollard/files/btag_ntuples/prod_ttbar/dy_msup50.root");
   TH1F* inputNumber3 = (TH1F*)file_dy2->Get("ttbarselectionproducer/hcheck_cutflow");
   n_dy2 = inputNumber3->GetBinContent(1);
   TChain* tree_dy2 = (TChain*) file_dy2->Get("btagana/ttree");
   CommPlotProducer4ttbar m_dy2(tree_dy2);


   TFile* file_st1 = new TFile("/opt/sbg/data/data3/cms/ccollard/files/btag_ntuples/prod_ttbar/single_t.root");
   TH1F* inputNumber4 = (TH1F*)file_st1->Get("ttbarselectionproducer/hcheck_cutflow");
   n_st1 = inputNumber4->GetBinContent(1);
   TChain* tree_st1 = (TChain*) file_st1->Get("btagana/ttree");
   CommPlotProducer4ttbar m_st1(tree_st1);

   TFile* file_st2 = new TFile("/opt/sbg/data/data3/cms/ccollard/files/btag_ntuples/prod_ttbar/single_tbar.root");
   TH1F* inputNumber5 = (TH1F*)file_st2->Get("ttbarselectionproducer/hcheck_cutflow");
   n_st2 = inputNumber5->GetBinContent(1);
   TChain* tree_st2 = (TChain*) file_st2->Get("btagana/ttree");
   CommPlotProducer4ttbar m_st2(tree_st2);


   // X-sections  (from AN2012_372_v2)

   float sf=0.978*0.986;
   float sf_dy=1.33;
   m_ttbar.SetNorm(227*sf*lumi/n_ttbar);   // 234 Theo, 227 mesure CMS TOP-12-007
   m_dy1.SetNorm(860.5*sf_dy*lumi/n_dy1);
   m_dy2.SetNorm(3532.8*sf_dy*lumi/n_dy2);
   m_st1.SetNorm(11.2*lumi/n_st1);
   m_st2.SetNorm(11.2*lumi/n_st2);

   // PU
   TString PUdataFile="../Production/lumiSum/MyDataPileupHistogram.root";
   SetPU2012_S10(PUdataFile);

   // Loop 
   m_data.Loop(0,3, 30, 470, "output_data");
   m_ttbar.Loop(1,3, 30, 470, "output_ttbar");
   m_dy1.Loop(2,3, 30, 470, "output_dy1");
   m_dy2.Loop(2,3, 30, 470, "output_dy2");
   m_st1.Loop(3,3, 30, 470, "output_st1");
   m_st2.Loop(3,3, 30, 470, "output_st2");

   // merge
}
