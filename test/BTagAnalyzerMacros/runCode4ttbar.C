{
   
   float luminosity = 2488; // /pb  see mbuttign@lxplus0058:~/BTagCommissioningPlots/CMSSW_7_6_3/src/RecoBTag/PerformanceMeasurements/test/ttbar and do brilcalc lumi -i lumiSummary_MuonEG_CD.json

  std::vector<TString > systlist;
  systlist.push_back(""                );
/*
  systlist.push_back("lept__plus"      );
  systlist.push_back("lept__minus"     );
  systlist.push_back("trig__plus"      );
  systlist.push_back("trig__minus"     );
  //systlist.push_back("PDF__plus"       );
  //systlist.push_back("PDF__minus"      );
  systlist.push_back("PU__plus"        );
  systlist.push_back("PU__minus"       );
  systlist.push_back("jes__plus"       );
  systlist.push_back("jes__minus"      );
  systlist.push_back("jer__plus"       );
  systlist.push_back("jer__minus"      );
  systlist.push_back("scale1__plus"    );
  systlist.push_back("scale1__minus"   );
  systlist.push_back("scale2__plus"    );
  systlist.push_back("scale2__minus"   );
*/

  //systlist.push_back("metuncls__plus"  );
  //systlist.push_back("metuncls__minus" );
  //systlist.push_back("toppt__plus"     );
  //systlist.push_back("toppt__minus"    );
  //systlist.push_back("btag__plus"      );
  //systlist.push_back("btag__minus"     );
  //systlist.push_back("mistag__plus"    );
  //systlist.push_back("mistag__minus"   );




   // initialisation
   TFile* file_data = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/Data13TeV_MuonEG_2015CD/MergedJetTree.root");
   TChain* tree_data = (TChain*) file_data->Get("btagana/ttree");

   double n_ttbar=0, wgtcounter_ttbar=0;
   double n_dy1=0, wgtcounter_dy1=0;
   double n_dy2=0, wgtcounter_dy2=0;
   double n_st1=0, wgtcounter_st1=0;
   double n_st2=0, wgtcounter_st2=0;
   double n_ww=0,  wgtcounter_ww=0;
   double n_wz=0,  wgtcounter_wz=0;
   double n_zz=0,  wgtcounter_zz=0;

   //TFile* file_ttbar = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_TTJets_Pedro/MergedJetTree.root");
   //TFile* file_ttbar = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_TTJets_madgraphMLM_pythia8/MergedJetTree.root");
   TFile* file_ttbar = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_TTJets_powheg_pythia8/MergedJetTree.root");
   //TFile* file_ttbar = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_TTJets_amcatnlo/MergedJetTree.root");
   TH1F* inputWeight = (TH1F*)file_ttbar->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_ttbar = inputWeight->GetBinContent(1);
   TChain* tree_ttbar = (TChain*) file_ttbar->Get("btagana/ttree");

   TFile* file_dy1 = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_DY10to50/MergedJetTree.root");
   TH1F* inputWeight2 = (TH1F*)file_dy1->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_dy1 = inputWeight2->GetBinContent(1);
   TChain* tree_dy1 = (TChain*) file_dy1->Get("btagana/ttree");

   TFile* file_dy2 = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_DY50toInf/MergedJetTree.root");
   TH1F* inputWeight3 = (TH1F*)file_dy2->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_dy2 = inputWeight3->GetBinContent(1);
   TChain* tree_dy2 = (TChain*) file_dy2->Get("btagana/ttree");

   TFile* file_st1 = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_SingleT_tW/MergedJetTree.root");
   TH1F* inputWeight4 = (TH1F*)file_st1->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_st1 = inputWeight4->GetBinContent(1);
   TChain* tree_st1 = (TChain*) file_st1->Get("btagana/ttree");

   TFile* file_st2 = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_SingleTbar_tW/MergedJetTree.root");
   TH1F* inputWeight5 = (TH1F*)file_st2->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_st2 = inputWeight5->GetBinContent(1);
   TChain* tree_st2 = (TChain*) file_st2->Get("btagana/ttree");

   TFile* file_ww = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_WWTo2L2Nu/MergedJetTree.root");
   TH1F* inputWeight6 = (TH1F*)file_ww->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_ww = inputWeight6->GetBinContent(1);
   TChain* tree_ww = (TChain*) file_ww->Get("btagana/ttree");

   TFile* file_wz = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_WZ/MergedJetTree.root");
   TH1F* inputWeight7 = (TH1F*)file_wz->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_wz = inputWeight7->GetBinContent(1);
   TChain* tree_wz = (TChain*) file_wz->Get("btagana/ttree");

   TFile* file_zz = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_ZZ/MergedJetTree.root");
   TH1F* inputWeight8 = (TH1F*)file_zz->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_zz = inputWeight8->GetBinContent(1);
   TChain* tree_zz = (TChain*) file_zz->Get("btagana/ttree");

   TFile* file_ttscaleup = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_TTJets_powheg_scaleup_pythia8/MergedJetTree.root");
   TH1F* inputWeight9 = (TH1F*)file_ttscaleup->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_ttscaleup = inputWeight9->GetBinContent(1);
   TChain* tree_ttscaleup = (TChain*) file_ttscaleup->Get("btagana/ttree");

   TFile* file_ttscaledown = new TFile("/opt/sbg/data/data2/cms/mbuttign/NTuplesBtagCommissioning_Jan16_v2/MC13TeV_TTJets_powheg_scaledown_pythia8/MergedJetTree.root");
   TH1F* inputWeight10 = (TH1F*)file_ttscaledown->Get("ttbarselectionproducer/wgtcounter");
   wgtcounter_ttscaledown = inputWeight10->GetBinContent(1);
   TChain* tree_ttscaledown = (TChain*) file_ttscaledown->Get("btagana/ttree");


   float sf_dy=1.3;


   CommPlotProducer4ttbar *m_data = new CommPlotProducer4ttbar(tree_data,1,1);
   m_data->Loop(0,3, 30, 500, "output_data_mueg", 0, "");


   // Loop on syst
   for (unsigned int syst = 0; syst < systlist.size(); syst++)
   { 

        if(systlist[syst] != "scale1__minus" && systlist[syst] != "scale1__plus" && systlist[syst] != "scale2__minus" && systlist[syst] != "scale2__plus")
        {

                CommPlotProducer4ttbar *m_ttbar = new CommPlotProducer4ttbar(tree_ttbar,1,1);
                m_ttbar->SetNorm(831.77*luminosity/wgtcounter_ttbar);  
                m_ttbar->Loop(1,3, 30, 500, "output_ttbar", inputWeight, systlist[syst]);

                CommPlotProducer4ttbar *m_dy1 = new CommPlotProducer4ttbar(tree_dy1,1,1);
                m_dy1->SetNorm(18610*sf_dy*luminosity/wgtcounter_dy1);
                m_dy1->Loop(2,3, 30, 500, "output_dy1", inputWeight2, systlist[syst]);

                CommPlotProducer4ttbar *m_dy2 = new CommPlotProducer4ttbar(tree_dy2,1,1);
                m_dy2->SetNorm(6025*sf_dy*luminosity/wgtcounter_dy2);
                m_dy2->Loop(2,3, 30, 500, "output_dy2", inputWeight3, systlist[syst]);

                CommPlotProducer4ttbar *m_st1 = new CommPlotProducer4ttbar(tree_st1,1,1);
                m_st1->SetNorm(35.85*luminosity/wgtcounter_st1);
                m_st1->Loop(3,3, 30, 500, "output_st1", inputWeight4, systlist[syst]);

                CommPlotProducer4ttbar *m_st2 = new CommPlotProducer4ttbar(tree_st2,1,1);
                m_st2->SetNorm(35.85*luminosity/wgtcounter_st2);
                m_st2->Loop(3,3, 30, 500, "output_st2", inputWeight5, systlist[syst]);

                CommPlotProducer4ttbar *m_ww = new CommPlotProducer4ttbar(tree_ww,1,1);
                m_ww->SetNorm(12.178*luminosity/wgtcounter_ww);
                m_ww->Loop(4,3, 30, 500, "output_ww", inputWeight6, systlist[syst]);

                CommPlotProducer4ttbar *m_wz = new CommPlotProducer4ttbar(tree_wz,1,1);
                m_wz->SetNorm(47.13*luminosity/wgtcounter_wz);
                m_wz->Loop(5,3, 30, 500, "output_wz", inputWeight7, systlist[syst]);

                CommPlotProducer4ttbar *m_zz = new CommPlotProducer4ttbar(tree_zz,1,1);
                m_zz->SetNorm(16.5*luminosity/wgtcounter_zz);
                m_zz->Loop(6,3, 30, 500, "output_zz", inputWeight8, systlist[syst]);


        }
        else if(systlist[syst] == "scale1__minus" || systlist[syst] == "scale1__plus")
        {
                CommPlotProducer4ttbar *m_ttbar = new CommPlotProducer4ttbar(tree_ttbar,1,1);
                m_ttbar->SetNorm(831.77*luminosity/wgtcounter_ttbar);  
                m_ttbar->Loop(1,3, 30, 500, "output_ttbar", inputWeight, systlist[syst]);
        }
        else if(systlist[syst] == "scale2__plus")
        {
                CommPlotProducer4ttbar *m_ttscaleup = new CommPlotProducer4ttbar(tree_ttscaleup,1,1);
                m_ttscaleup->SetNorm(831.77*luminosity/wgtcounter_ttscaleup);
                m_ttscaleup->Loop(1,3, 30, 500, "output_ttscaleup", inputWeight9, systlist[syst]);
        }
        else if(systlist[syst] == "scale2__minus")
        {
                CommPlotProducer4ttbar *m_ttscaledown = new CommPlotProducer4ttbar(tree_ttscaledown,1,1);
                m_ttscaledown->SetNorm(831.77*luminosity/wgtcounter_ttscaledown);
                m_ttscaledown->Loop(1,3, 30, 500, "output_ttscaledown", inputWeight10, systlist[syst]);
        }       

   }       

}
