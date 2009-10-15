{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L JetsAna.C+g") ;
  
 if (gROOT->GetClass("Jets")==0) return;
 
 TChain c("Jets");

////////////////////////////////////////////////////////////////////////////////

 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt30_7TeV_SD_Jet50U_1.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_1.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_2.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_3.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_4.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_5.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_6.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_7.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_8.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_9.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_10.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_11.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_12.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_14.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_13.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_14.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_15.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_16.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_17.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_18.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_19.root");
 c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_7TeV_SD_Jet50U_20.root");

//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt30_SD_Jet50U_1.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_1.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_2.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_3.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_4.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_5.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_6.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_7.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_8.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_9.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_10.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_11.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_12.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_13.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/Ntup_QCD_Pt80_SD_Jet50U_14.root");

//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_1.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_2.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_3.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_4.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_5.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_6.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_7.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_8.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_9.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_10.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_11.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt80_7TeV_12.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_1.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_2.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_3.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_4.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_5.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_6.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_7.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_8.root");
//  c.Add("../NTUPLES/CMSSW_3_1_4/NtupTkHist_QCDpt170_7TeV_9.root");

////////////////////////////////////////////////////////////////////////////////

 Jets* t = new Jets(Jets);

//  t->Loop(0,0.50,30.,999.,0.,2.4,1000.,"output.root");

 t->Loop(0,.230,30.,999.,0.0,0.7,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_JPL.root");
 t->Loop(0,.230,30.,999.,0.7,1.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_JPL.root");
 t->Loop(0,.230,30.,999.,1.4,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_JPL.root");
 t->Loop(0,.230,30.,999.,0.0,0.7,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_JPL_veto.root");
 t->Loop(0,.230,30.,999.,0.7,1.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_JPL_veto.root");
 t->Loop(0,.230,30.,999.,1.4,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_JPL_veto.root");

 t->Loop(0,.495,30.,999.,0.0,0.7,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_JPM.root");
 t->Loop(0,.495,30.,999.,0.7,1.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_JPM.root");
 t->Loop(0,.495,30.,999.,1.4,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_JPM.root");
 t->Loop(0,.495,30.,999.,0.0,0.7,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_JPM_veto.root");
 t->Loop(0,.495,30.,999.,0.7,1.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_JPM_veto.root");
 t->Loop(0,.495,30.,999.,1.4,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_JPM_veto.root");

 t->Loop(0,.700,30.,999.,0.0,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_JPT.root");
 t->Loop(0,.700,30.,999.,0.0,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_JPT_veto.root");

 t->Loop(2,1.90,30.,999.,0.0,0.7,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_TCHEL.root");
 t->Loop(2,1.90,30.,999.,0.7,1.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_TCHEL.root");
 t->Loop(2,1.90,30.,999.,1.4,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_TCHEL.root");
 t->Loop(2,1.90,30.,999.,0.0,0.7,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_TCHEL_veto.root");
 t->Loop(2,1.90,30.,999.,0.7,1.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_TCHEL_veto.root");
 t->Loop(2,1.90,30.,999.,1.4,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_TCHEL_veto.root");

 t->Loop(2,3.99,30.,999.,0.0,0.7,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_TCHEM.root");
 t->Loop(2,3.99,30.,999.,0.7,1.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_TCHEM.root");
 t->Loop(2,3.99,30.,999.,1.4,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_TCHEM.root");
 t->Loop(2,3.99,30.,999.,0.0,0.7,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_TCHEM_veto.root");
 t->Loop(2,3.99,30.,999.,0.7,1.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_TCHEM_veto.root");
 t->Loop(2,3.99,30.,999.,1.4,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_TCHEM_veto.root");

 t->Loop(3,2.17,30.,999.,0.0,0.7,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_TCHPM.root");
 t->Loop(3,2.17,30.,999.,0.7,1.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_TCHPM.root");
 t->Loop(3,2.17,30.,999.,1.4,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_TCHPM.root");
 t->Loop(3,2.17,30.,999.,0.0,0.7,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_TCHPM_veto.root");
 t->Loop(3,2.17,30.,999.,0.7,1.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_TCHPM_veto.root");
 t->Loop(3,2.17,30.,999.,1.4,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_TCHPM_veto.root");

 t->Loop(3,4.31,30.,999.,0.0,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_TCHPT.root");
 t->Loop(3,4.31,30.,999.,0.0,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_TCHPT_veto.root");

 t->Loop(4,2.02,30.,999.,0.0,0.7,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_SSVM.root");
 t->Loop(4,2.02,30.,999.,0.7,1.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_SSVM.root");
 t->Loop(4,2.02,30.,999.,1.4,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_SSVM.root");
 t->Loop(4,2.02,30.,999.,0.0,0.7,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_etaLT07_SSVM_veto.root");
 t->Loop(4,2.02,30.,999.,0.7,1.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_07eta14_SSVM_veto.root");
 t->Loop(4,2.02,30.,999.,1.4,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_14eta24_SSVM_veto.root");

 t->Loop(4,3.40,30.,999.,0.0,2.4,1000.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_SSVT.root");
 t->Loop(4,3.40,30.,999.,0.0,2.4,20.,"output/QCD_Pt30-80_7TeV_SD_Jet50U_SSVT_veto.root");

 gSystem->Exec("kill -9 "+TString(Form("%d",gSystem->GetPid())));

////////////////////////////////////////////////////////////////////////////////
}
