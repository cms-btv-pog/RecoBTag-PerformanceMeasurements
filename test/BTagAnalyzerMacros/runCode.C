{
    gROOT->ProcessLine(".L CommPlotProducer.C++");
   
    // Declare the root files on which you want to run  (MC & Data) :
    TChain *superTree = new TChain("btagana/ttree");

    int config_ch = 2;
    cout << " which config? Phys14=1, CSA14=2, 8TeV=3 " << endl;
    cin >> config_ch ;
    cout << " config choisie " << config_ch << endl;
    // MC
    if (config_ch==1) {
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/Phys14_wtrk/QCD_Pt-30to50_Tune4C_13TeV_pythia8/0000/JetTree_phys14_*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/Phys14_wtrk/QCD_Pt-50to80_Tune4C_13TeV_pythia8/0000/JetTree_phys14_*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/Phys14_wtrk/QCD_Pt-80to120_Tune4C_13TeV_pythia8/0000/JetTree_phys14_*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/Phys14_wtrk/QCD_Pt-120to170_Tune4C_13TeV_pythia8/0000/JetTree_phys14_*.root");

    }
    else if (config_ch==2) {
      // CSA14
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-30to50_Tune4C_13TeV_pythia8/0000/JetTree_csa14*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-50to80_Tune4C_13TeV_pythia8/0000/JetTree_csa14*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-80to120_Tune4C_13TeV_pythia8/0000/JetTree_csa14*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-120to170_Tune4C_13TeV_pythia8/0000/JetTree_csa14*.root");
    }
    else if (config_ch==3) {
      // 8 TeV
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/QCD_Pt-30to50_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/QCD_Pt-30to50_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/JetTree_mc_*root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/QCD_Pt-50to80_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/JetTree_mc_*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/QCD_Pt-80to120_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/JetTree_mc_*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/QCD_Pt-120to170_Inclusive_8TeV_pythia6_Summer12_DR53X-PU_S10_AODSIM/JetTree_mc_*.root");

    }
    else if (config_ch==4) {
      // CSA14 Pythia6
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-30to50_TuneZ2star_13TeV_pythia6/0000/JetTree_csa14*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-50to80_TuneZ2star_13TeV_pythia6/0000/JetTree_csa14*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-80to120_TuneZ2star_13TeV_pythia6/0000/JetTree_csa14*.root");
      superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/CSA14/newProd/QCD_Pt-120to170_TuneZ2star_13TeV_pythia6/0000/JetTree_csa14*.root");

    }
    else if (config_ch==5) {
      // 8 Tev Pythia8
      superTree->Add("/opt/sbg/data/data2/cms/ccollard/BTAG/prod_pythia8/pythia8_hadd/JetTree_mc_qcd30.root");
      superTree->Add("/opt/sbg/data/data2/cms/ccollard/BTAG/prod_pythia8/pythia8_hadd/JetTree_mc_qcd50.root");
      superTree->Add("/opt/sbg/data/data2/cms/ccollard/BTAG/prod_pythia8/pythia8_hadd/JetTree_mc_qcd80.root");
      superTree->Add("/opt/sbg/data/data2/cms/ccollard/BTAG/prod_pythia8/pythia8_hadd/JetTree_mc_qcd120.root");
    }


      CommPlotProducer m(superTree);

    // Provide MC information. use SetXS(TString generator, bool qcdtype, int TeV) ;
    // with generator = pythia or herwig, qcdtype =0 for inclusive or 1 for MuEnriched, and TeV = 8 or 13.
    if (config_ch!=3 && config_ch<5) {
        m.SetInfo("pythia",0,13);
    }
    else {
        m.SetInfo("pythia",0,8);
    }
 
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
    double   n30    = 0.; 
    double   n50    = 0.;
    double   n80    = 0.;
    double   n120   = 0.; 
    double   n170   = 0.;
    double   n300   = 0.;
    double   n470   = 0.; 
    double   n600   = 0.;  
    double   n800   = 0.;  
    double   n1000  = 0.;

    if (config_ch==1) { 
        n30   = 498559; 
        n50   = 517391; 
        n80   = 924531; 
        n120   = 1015701; 
    }
    else if (config_ch==2) { 
        n30   = 861636; 
        n50   = 831797; 
        n80   = 866767; 
        n120   = 746130; 
    }
    else if (config_ch==3) { 
        n30   = 5998027; 
        n50   = 5994786; 
        n80   = 5986846; 
        n120   = 5970933; 
    }
    else if (config_ch==4) {
        n30   = 547116; 
        n50   = 588659;
        n80   = 600000;
        n120   = 596768;
    }
    else if (config_ch==5) {
        n30   = 1000080;
        n50    = 1000026;
        n80   =  1000054;
        n120  =   800064;
    }


    m.Fill_nevent(0.,n20,n30,n50,n80,n120,n170,n300,n470,n600,n800,n1000);
   

    // Set.XS() will automatically put the cross sections depending on the choice make above with SetInfo().
    // if nothing was defined, the default is to use inclusive pythia x-sections for 8 TeV.
    m.SetXS();
 
    // Compute the Total x-section of all the samples
    m.SetSumXS();  
   
    // Here put the name of the PU data root file for the PU reweighting
    //TString PUdataFile="PileUp_jet_abcd.root";

    // For the PU reweighting, you can use the 2012 MC distribution with PU sceanrio S7 (SetPU2012_S7(PUdataFile)),
    // with PU scenario S10 (SetPU2012_S10(PUdataFile)), or with your own PU vector. Be carefull this vector has to be
    // of size 60.
   
    //m.SetPU2012_S10(PUdataFile);
   
    //The loop on tthe events takes 4 inputs: the trigger, the ptmin of the jets, the ptmax of the jets, and the output
    //file name 
   
    TString name_root;
    if (config_ch==1) {  name_root = "phys14"; }
    else if (config_ch==2) {  name_root = "csa14"; }
    else if (config_ch==3) {  name_root = "8tev_new"; }
    else if (config_ch==4) {  name_root = "csa14_pythia6"; }
    else if (config_ch==5) {  name_root = "8tev_pythia8"; }
    m.Loop("jet", 20, 30, 200,   name_root);
   
}
