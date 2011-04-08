{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L JetsAna.C+g") ;
  
 if (gROOT->GetClass("JetTree")==0) return;
 
 TChain c("mistag/ttree");

////////////////////////////////////////////////////////////////////////////////

c.Add("../NTUPLES/CMSSW_4_1_4/JetTree.root");

// c.Add("../NTUPLES/CMSSW_4_1_2_patch1/Data/Jet_Run2011A-PromptReco/JetTree_160404-161312.root");

// c.Add("../NTUPLES/CMSSW_4_1_2_patch1/MC/QCD_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1/JetTree_Pt_15to30.root"); // Jet30
// c.Add("../NTUPLES/CMSSW_4_1_2_patch1/MC/QCD_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1/JetTree_Pt_30to50.root"); // Jet60
// c.Add("../NTUPLES/CMSSW_4_1_2_patch1/MC/QCD_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1/JetTree_Pt_50to80.root"); // Jet80
// c.Add("../NTUPLES/CMSSW_4_1_2_patch1/MC/QCD_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1/JetTree_Pt_80to120.root");  // Jet110
// c.Add("../NTUPLES/CMSSW_4_1_2_patch1/MC/QCD_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1/JetTree_Pt_120to170.root"); // Jet150 + 190
// c.Add("../NTUPLES/CMSSW_4_1_2_patch1/MC/QCD_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1/JetTree_Pt_170to300.root");
// // c.Add("../NTUPLES/CMSSW_4_1_2_patch1/MC/QCD_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1/JetTree_Pt_300to470.root");

////////////////////////////////////////////////////////////////////////////////

 JetTree* t = new JetTree(&c);

 t->Loop(2,1.70,20.,999.,0.,2.4,0.,0,"output.root");

// //  t->Loop(0,.275,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_JPL.root");
//  t->Loop(0,.545,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_JPM.root");
// //  t->Loop(0,.790,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_JPT.root");
// //  t->Loop(1,1.33,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_JBPL.root");
//  t->Loop(1,2.55,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_JBPM.root");
// //  t->Loop(1,3.74,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_JBPT.root");
// //  t->Loop(2,1.70,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_TCHEL.root");
//  t->Loop(2,3.30,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_TCHEM.root");
//  t->Loop(3,1.93,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_TCHPM.root");
// //  t->Loop(3,3.41,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_TCHPT.root");
//  t->Loop(4,1.74,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_SSVHEM.root");
//  t->Loop(5,2.00,20.,999.,0.,2.4,0.,30,"output/QCD_Jet30_SSVHPT.root");

// //  t->Loop(0,.275,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_JPL.root");
//  t->Loop(0,.545,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_JPM.root");
// //  t->Loop(0,.790,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_JPT.root");
// //  t->Loop(1,1.33,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_JBPL.root");
//  t->Loop(1,2.55,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_JBPM.root");
// //  t->Loop(1,3.74,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_JBPT.root");
// //  t->Loop(2,1.70,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_TCHEL.root");
//  t->Loop(2,3.30,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_TCHEM.root");
//  t->Loop(3,1.93,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_TCHPM.root");
// //  t->Loop(3,3.41,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_TCHPT.root");
//  t->Loop(4,1.74,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_SSVHEM.root");
//  t->Loop(5,2.00,20.,999.,0.,2.4,0.,60,"output/QCD_Jet60_SSVHPT.root");
// 
// //  t->Loop(0,.275,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_JPL.root");
//  t->Loop(0,.545,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_JPM.root");
// //  t->Loop(0,.790,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_JPT.root");
// //  t->Loop(1,1.33,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_JBPL.root");
//  t->Loop(1,2.55,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_JBPM.root");
// //  t->Loop(1,3.74,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_JBPT.root");
// //  t->Loop(2,1.70,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_TCHEL.root");
//  t->Loop(2,3.30,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_TCHEM.root");
//  t->Loop(3,1.93,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_TCHPM.root");
// //  t->Loop(3,3.41,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_TCHPT.root");
//  t->Loop(4,1.74,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_SSVHEM.root");
//  t->Loop(5,2.00,20.,999.,0.,2.4,0.,80,"output/QCD_Jet80_SSVHPT.root");
// 
// //  t->Loop(0,.275,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_JPL.root");
//  t->Loop(0,.545,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_JPM.root");
// //  t->Loop(0,.790,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_JPT.root");
// //  t->Loop(1,1.33,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_JBPL.root");
//  t->Loop(1,2.55,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_JBPM.root");
// //  t->Loop(1,3.74,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_JBPT.root");
// //  t->Loop(2,1.70,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_TCHEL.root");
//  t->Loop(2,3.30,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_TCHEM.root");
//  t->Loop(3,1.93,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_TCHPM.root");
// //  t->Loop(3,3.41,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_TCHPT.root");
//  t->Loop(4,1.74,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_SSVHEM.root");
//  t->Loop(5,2.00,20.,999.,0.,2.4,0.,110,"output/QCD_Jet110_SSVHPT.root");
// 
// //  t->Loop(0,.275,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_JPL.root");
//  t->Loop(0,.545,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_JPM.root");
// //  t->Loop(0,.790,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_JPT.root");
// //  t->Loop(1,1.33,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_JBPL.root");
//  t->Loop(1,2.55,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_JBPM.root");
// //  t->Loop(1,3.74,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_JBPT.root");
// //  t->Loop(2,1.70,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_TCHEL.root");
//  t->Loop(2,3.30,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_TCHEM.root");
//  t->Loop(3,1.93,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_TCHPM.root");
// //  t->Loop(3,3.41,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_TCHPT.root");
//  t->Loop(4,1.74,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_SSVHEM.root");
//  t->Loop(5,2.00,20.,999.,0.,2.4,0.,150,"output/QCD_Jet150_SSVHPT.root");
// 
// //  t->Loop(0,.275,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_JPL.root");
//  t->Loop(0,.545,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_JPM.root");
// //  t->Loop(0,.790,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_JPT.root");
// //  t->Loop(1,1.33,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_JBPL.root");
//  t->Loop(1,2.55,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_JBPM.root");
// //  t->Loop(1,3.74,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_JBPT.root");
// //  t->Loop(2,1.70,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_TCHEL.root");
//  t->Loop(2,3.30,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_TCHEM.root");
//  t->Loop(3,1.93,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_TCHPM.root");
// //  t->Loop(3,3.41,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_TCHPT.root");
//  t->Loop(4,1.74,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_SSVHEM.root");
//  t->Loop(5,2.00,20.,999.,0.,2.4,0.,190,"output/QCD_Jet190_SSVHPT.root");

//  t->Loop(2,1.70,20.,999.,0.,2.4,0.,0,"output/QCD_TCHEL.root");

////////////////////////////////////////////////////////////////////////////////

//  gSystem->Exec("kill -9 "+TString(Form("%d",gSystem->GetPid())));
}
