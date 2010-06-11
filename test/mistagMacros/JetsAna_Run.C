{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L JetsAna.C+g") ;
  
 if (gROOT->GetClass("JetTree")==0) return;
 
 TChain c("mistag/ttree");

////////////////////////////////////////////////////////////////////////////////

c.Add("~/NTUPLES/CMSSW_3_6_1_patch3/Data/JetMETTau_Run2010A-PromptReco-v2/JetTree_All.root");

////////////////////////////////////////////////////////////////////////////////

 JetTree* t = new JetTree(&c);

 t->Loop(5,2.00,30.,999.,0.,2.4,1000.,"Data-prompt-Jet15U_ptGT30_SSVHPT.root");
//  t->Loop(2,1.25,30.,999.,0.,2.4,1000.,"output.root");

//  t->Loop(1,3.00,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_IP1.root");
//  t->Loop(2,1.70,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_TCHEL.root");
//  t->Loop(2,3.30,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_TCHEM.root");
//  t->Loop(3,1.93,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_TCHPM.root");
//  t->Loop(3,3.41,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_TCHPT.root");
//  t->Loop(4,1.74,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_SSVHEM.root");
//  t->Loop(4,3.05,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_SSVHET.root");
//  t->Loop(5,2.00,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_SSVHPT.root");
//  t->Loop(0,.215,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_JPL.root");
//  t->Loop(0,.459,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_JPM.root");
//  t->Loop(0,.669,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_JPT.root");
//  t->Loop(5,0.38,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_CSVL.root");
//  t->Loop(5,0.75,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_CSVM.root");
//  t->Loop(5,.921,0.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_CSVT.root");

//  t->Loop(1,3.00,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_IP1.root");
//  t->Loop(2,1.70,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_TCHEL.root");
//  t->Loop(2,3.30,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_TCHEM.root");
//  t->Loop(3,1.93,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_TCHPM.root");
//  t->Loop(3,3.41,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_TCHPT.root");
//  t->Loop(4,1.74,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_SSVHEM.root");
//  t->Loop(4,3.05,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_SSVHET.root");
//  t->Loop(5,2.00,30.,999.,0.,2.4,1000.,"output/Data-prompt-Jet15U_ptGT30_SSVHPT.root");

////////////////////////////////////////////////////////////////////////////////

 gSystem->Exec("kill -9 "+TString(Form("%d",gSystem->GetPid())));
}
