{
   
    // Declare the root files on which you want to run  (MC & Data) :
    TChain *superTree = new TChain("btagana/ttree");

    int config_ch = 2;
    cout << " which config? from 1 to 6 for QCD and from 7 to 9 for TTbar " << endl;
    cin >> config_ch ;
    cout << " config choisie " << config_ch << endl;
    // MC
    if (config_ch==1) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/QCD/QCD50to80_8TeVPythia6_NoPU.root");
    else if (config_ch==2) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/QCD/QCD50to80_8TeVPythia6_PU.root");
    else if (config_ch==3) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/QCD/QCD50to80_8TeVPythia8_NoPU.root");
    else if (config_ch==4) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/QCD/QCD50to80_8TeVPythia8_PU.root");
    else if (config_ch==5) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/QCD/QCD50to80_13TeVPythia8_NoPU.root");
    else if (config_ch==6) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/QCD/QCD50to80_13TeVPythia8_PU.root");
    else if (config_ch==7) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/TTbar/TTbar_8TeVPythia6_PU.root");
    else if (config_ch==8) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/TTbar/TTbar_8TeVPythia8_NoPU.root");
    else if (config_ch==9) superTree->Add("/opt/sbg/data/data5/cms/ccollard/BTag/BTagAnalyzerNtuples/RunIWinter15DR/TTbar/TTbar_13TeVPythia8_NoPU.root");


    CommPlotProducer m(superTree);


    // Setting up of the information
    // with generator = pythia or herwig, qcdtype =0 for inclusive or 1 for MuEnriched, and TeV = 8 or 13.
    if (config_ch<5) m.SetInfo("pythia",0,8);
    else if (config_ch==5 || config_ch==6)  m.SetInfo("pythia",0,13);
    else m.SetInfo("pythia",0,-1);
 
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

   
    if (config_ch==1)      m.Fill_nevent(0.,0.,0., 99360.,0.,0.,0.,0.,0.,0.,0.,0.);
    else if (config_ch==2) m.Fill_nevent(0.,0.,0., 99360.,0.,0.,0.,0.,0.,0.,0.,0.);
    else if (config_ch==3) m.Fill_nevent(0.,0.,0.,100000.,0.,0.,0.,0.,0.,0.,0.,0.);
    else if (config_ch==4) m.Fill_nevent(0.,0.,0.,100000.,0.,0.,0.,0.,0.,0.,0.,0.);
    else if (config_ch==5) m.Fill_nevent(0.,0.,0.,100000.,0.,0.,0.,0.,0.,0.,0.,0.);
    else if (config_ch==6) m.Fill_nevent(0.,0.,0.,100000.,0.,0.,0.,0.,0.,0.,0.,0.);

    if (config_ch<7)  {
      m.SetXS();     // Assign the correct x-sections to QCD pthat bins, depending on SetInfo(), default = use inclusive pythia x-sections for 8 TeV.
      m.SetSumXS();  // Compute the Total x-section of all the samples
    }
   
    // no PU info
   
    TString name_root;
    if (config_ch==1)      {  name_root = "QCD_8P6NoPU"; }
    else if (config_ch==2) {  name_root = "QCD_8P6wPU"; }
    else if (config_ch==3) {  name_root = "QCD_8P8NoPU"; }
    else if (config_ch==4) {  name_root = "QCD_8P8wPU"; }
    else if (config_ch==5) {  name_root = "QCD_13P8NoPU"; }
    else if (config_ch==6) {  name_root = "QCD_13P8wPU"; }
    else if (config_ch==7) {  name_root = "TT_8P6wPU"; }
    else if (config_ch==8) {  name_root = "TT_8P8NoPU"; }
    else if (config_ch==9) {  name_root = "TT_13P8NoPU"; }

    m.Loop("jet", 20, 30, 80,   name_root);
   
}
