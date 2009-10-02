{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L JetsAna.C+g") ;
  
 if (gROOT->GetClass("Jets")==0) return;
 
 TChain c("Jets");

////////////////////////////////////////////////////////////////////////////////

 c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_skimHLT_Jet50U_1.root");

//  c.Add("~/NTUPLES/CMSSW_3_1_2/NtupTkHist_QCDpt300_preprod_1.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/NtupTkHist_QCDpt300_preprod_2.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/NtupTkHist_QCDpt300_preprod_3.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/NtupTkHist_QCDpt300_preprod_4.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/NtupTkHist_QCDpt300_preprod_5.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/NtupTkHist_QCDpt300_preprod_6.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/NtupTkHist_QCDpt300_preprod_7.root");

//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_1.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_2.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_3.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_4.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_5.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_6.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_7.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_8.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_9.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_10.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_11.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_12.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_13.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_14.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_15.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_16.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_17.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_18.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_19.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_20.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_21.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_22.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_23.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_24.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_25.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_26.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_27.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_28.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_29.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_30.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_31.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_32.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_33.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_34.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_35.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_36.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_37.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_38.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt30_39.root");

//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_2.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_3.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_4.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_5.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_7.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_8.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_9.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_10.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_11.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_13.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_14.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_15.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_16.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_17.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_18.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_19.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_20.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_21.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_22.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_24.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_25.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_26.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_27.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_28.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_29.root");
//  c.Add("~/NTUPLES/CMSSW_3_1_2/Ntup_QCDpt80_30.root");

////////////////////////////////////////////////////////////////////////////////

 Jets* t = new Jets(Jets);

//  t->Loop(0,0.50,30.,999.,0.,2.4,1000.,"output.root");

 t->Loop(0,.230,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_JPL.root");
 t->Loop(0,.495,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_JPM.root");
 t->Loop(0,.700,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_JPT.root");
 t->Loop(2,1.90,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_TCHEL.root");
 t->Loop(2,3.99,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_TCHEM.root");
 t->Loop(3,2.17,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_TCHPM.root");
 t->Loop(3,4.31,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_TCHPT.root");
 t->Loop(4,2.02,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_SSVM.root");
 t->Loop(4,3.40,30.,999.,0.,0.7,1000.,"output/QCDpt30_SDJet50U_etaLT07_SSVT.root");

 gSystem->Exec("kill -9 "+TString(Form("%d",gSystem->GetPid())));

////////////////////////////////////////////////////////////////////////////////
}
