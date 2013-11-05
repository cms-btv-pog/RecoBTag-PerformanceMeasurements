{
   gROOT->ProcessLine(".L CommPlotProducer.C++");
   
   // Declare the root files on which you want to run  (MC & Data) :
   TChain *superTree = new TChain("btagana/ttree");
//   superTree->Add("/opt/sbg/cms/ui2_data1/blochd/NTUPLES/CMSSW_5_3_11_patch1/MC/QCD_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A_V02-02-07/JetTree_subjets_Pt-50to80_403on415.root");

   // Data

   superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/Data/Jet_Run2012A-22Jan2013/*root");
   superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/Data/JetMon_Run2012B-22Jan2013/*root");
   superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/Data/JetMon_Run2012C-22Jan2013/*root");
   superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/Data/JetMon_Run2012D-22Jan2013/*root");

   // MC
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-30to50_TP.root");
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-50to80_TP.root");
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-80to120_TP.root");
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-120to170_TP.root");
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-170to300_TP.root");
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-300to470_TP.root");
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-470to600_TP.root");
  superTree->Add("/opt/sbg/cms/ui2_data2/blochd/NTUPLES/CMSSW_5_3_11_V02-03-02/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A/JetTree_mc_subjets_Pt-600to800_TP.root");

   // JetProbaTree : Oui
   // NewAlgoTree : Non
   CommPlotProducer m(superTree,0,1);

   // Provide MC information. use SetXS(TString generator, bool qcdtype, int TeV) ;
   // with generator = pythia or herwig, qcdtype =0 for inclusive or 1 for MuEnriched, and TeV =7 or 8.
   m.SetInfo("pythia",0,8);
 
   // Here put the number of events for each samples you are running on.
   // BE CAREFULL !!!!!
   // You should enter 12 numbers and following the increasing order in pthat of the MC samples.
   // n15_20, n20_30, n30_50, n50_80,  n80_120,  n120_170, n170_300,  n300_470,  n470_600, n600_800, n800_1000, n1000 ;
   //
   //  -->  m.Fill_nevent(n15_20,n20_30,n30_50,n50_80,n80_120,n120_170,n170_300,n300_470,n470_600,n600_800,n800_1000,n1000);
   //       for muenriched pythia and herwig.
   //  or
   //  -->  m.Fill_nevent(    0.,n15_30,n30_50,n50_80,n80_120,n120_170,n170_300,n300_470,n470_600,n600_800,n800_1000,n1000);
   //       for inclusive pythia and herwig.
   //
   // If you want to run only on the 50-80 and 80-120 samples for pythia 8 TeV, the correct syntaxe is
   //   m.Fill_nevent(0,0,0,n50_80,n80_120,0,0,0,0,0,0,0);
   //
   // Remark for 7TeV QCD MC:
   // for  inclusive QCD pythia at 7TeV, provide : 
   //  -->  m.Fill_nevent(    0.,n15_30,n30_50,n50_80,n80_120,n120_170,n170_300,n300,0.,0.,0.,0.);
   // for  MuEnriched QCD pythia at 7TeV, provide : 
   //  -->  m.Fill_nevent(    0.,n15_30,n30_50,n50_80,n80_120,n120_150,n150,0.,0.,0.,0.,0.);
   //
   // if you don't know how many events you have in your datasets
   // m.Counter(); --> TAKES 4 HOURS ON DATA+MC :(

    double   n20    = 0.;
    double   n30    = 408047; 
    double   n50    = 380139;
    double   n80   =  841606;
    double   n120  =  878681;
    double   n170  =  293229;
    double   n300  =  312488;
    double   n470  =  319603;
    double   n600  =  265460;
    double   n800 =  0.;
    double   n1000     = 0;
    m.Fill_nevent(0.,n20,n30,n50,n80,n120,n170,n300,n470,n600,n800, 0.);

   // Set.XS() will automatically put the cross sections depending on the choice make above with SetInfo().
   // if nothing was defined, the default is to use inclusive pythia x-sections for 8 TeV.
   m.SetXS();
 
   // Compute the Total x-section of all the samples
   m.SetSumXS();  
   
   // Here put the name of the PU data root file for the PU reweighting
   TString PUdataFile="PileUp_jet_abcd.root";

   // For the PU reweighting, you can use the 2012 MC distribution with PU sceanrio S7 (SetPU2012_S7(PUdataFile)),
   // with PU scenario S10 (SetPU2012_S10(PUdataFile)), or with your own PU vector. Be carefull this vector has to be
   // of size 60.
   
   m.SetPU2012_S10(PUdataFile);
   
   //The loop on tthe events takes 4 inputs: the trigger, the ptmin of the jets, the ptmax of the jets, and the output
   //file name 
   
   m.Loop("jet",80, 100, 500, "output_test80");
//   m.Loop("btag",20, 45, 300, "output_btag_newalgo1");


   
}
