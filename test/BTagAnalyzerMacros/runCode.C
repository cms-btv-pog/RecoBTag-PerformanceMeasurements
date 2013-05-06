{
   gROOT->ProcessLine(".L CommPlotProducer.C++");
   
   // Declare the root files on which you want to run  (MC & Data) :
   TChain *superTree = new TChain("btagana/ttree");
   // Data
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/Data_Winter13/Jet/Run2012A-v1/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/Data_Winter13/Jet/Run2012B-v1/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/Data_Winter13/Jet/Run2012C-v1/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/Data_Winter13/Jet/Run2012D-v1/TrackTree_*.root");
   // MC
   //superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-15to30_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-30to50_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-50to80_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-80to120_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data2/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-120to170_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data2/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-170to300_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-300to470_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data1/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-470to600_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data2/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-600to800_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");
   superTree->Add("/opt/sbg/cms/ui8_data2/ccollard/BTAG/MC_Dan_2013April/PythiaQCD/QCD_Pt-800to1000_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/TrackTree_*.root");

   // JetProbaTree : Oui
   // NewAlgoTree : Non
   CommPlotProducer m(superTree,1,0);

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

    double   n15    = 6.97142e+06;
    double   n20    = 2.93092e+06;
    double   n15to30 = n15+n20;
    double   n30    = 5.96571e+06;
    double   n50    = 5.98172e+06;
    double   n80   = 5.89208e+06;
    double   n120  = 5.91732e+06;
    double   n170  = 5.81438e+06;
    double   n300  = 5.89309e+06;
    double   n470  = 3.93778e+06;
    double   n600  = 3.98544e+06;
    double   n800 = 3.98713e+06;
    double   n1000     = 0;
    m.Fill_nevent(0.,n15to30,n30,n50,n80,n120,n170,n300,n470,n600,n800,n1000);
   // m.Fill_nevent(0.,0.,n30,n50,n80,n120,n170,n300,n470,0.,0., 0.);

   

   // Set.XS() will automatically put the cross sections depending on the choice make above with SetInfo().
   // if nothing was defined, the default is to use inclusive pythia x-sections for 8 TeV.
   m.SetXS();
 
   // Compute the Total x-section of all the samples
   m.SetSumXS();  
   
   // Here put the name of the PU data root file for the PU reweighting
   TString PUdataFile="../Production_data/lumiSum/PileUp_jet_abcd.root";

   // For the PU reweighting, you can use the 2012 MC distribution with PU sceanrio S7 (SetPU2012_S7(PUdataFile)),
   // with PU scenario S10 (SetPU2012_S10(PUdataFile)), or with your own PU vector. Be carefull this vector has to be
   // of size 60.
   
   m.SetPU2012_S10(PUdataFile);
   
   //The loop on tthe events takes 4 inputs: the trigger, the ptmin of the jets, the ptmax of the jets, and the output
   //file name 
   
   m.Loop("jet",40, 60, 500, "output_jet");

   
}
