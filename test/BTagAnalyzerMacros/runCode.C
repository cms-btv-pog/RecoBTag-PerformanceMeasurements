{
   gROOT->ProcessLine(".L CommPlotProducer.C++");
   
   // Declare the root files on which you want to run  (MC & Data) :
   TChain *superTree = new TChain("btagana/ttree");
   superTree->Add("/opt/sbg/data/data3/cms/ccollard/test_camille/CMSSW_5_3_2_patch4/src/RecoBTag/PerformanceMeasurements/test/TrackTree_*.root");
   CommPlotProducer m(superTree);

   // Here put the number of events for each samples you are running on.
   // BE CAREFULL !!!!!
   // You should enter 10 numbers and following the increasing order in pthat of the MC samples.
   // Depending you MC generator, use : 
   //    pythia 8TeV:             n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, 0;
   //    pythia 8 TeV Muenriched: n30-50, n50-80, n80-120, n120-170, n170-300, n300-470,  n470-600,  n600-800, 0       , 0;
   //    pythia 7TeV:             n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000;
   //    pythia 7TeV Mu-enriched: n15-20, n20-30, n30-50,  n50-80,   n80-120,  n120-150,  n150-plus, 0       , 0       , 0;
   //    herwig 8 TeV             n30-50, n50-80, n80-120, n120-170, n170-300, n300-plus, 0        , 0       , 0       , 0;
   // For example if you are using Pythia MC at 8 TeV with pthat bin 15-30, 30-50, 50-80,80-120,120-170,170-300,300-470,470-600 and 600-800 
   // the correct syntaxe will be
   //    m.Fill_nevent(n15-30,n30-50,n50-80,n80-120,n120-170,n170-300,n300-470,n470-600,n600-800,0);
   // If you want to run only on the 50-80 and 80-120 samples for pythia 8 TeV, the correct syntaxe is
   //   m.Fill_nevent(0,0,n50-80,n80-120,0,0,0,0,0,0);
   
   
   m.Fill_nevent(0,0,5989860.0,4356690.0,4031390.0,4031390.0,5804118.0,5972200.0,0,0);  
   // if you don't know how many events you have in your datasets
   //m.Counter();
   
   // Set.XS() will automatically put the cross sections for pythia  TeV. If you want to use a different generator at
   // a different energy, with muons samples, use SetXS(TString generator, bool MuEnriched, int TeV) ;
   // with generator = pythia or herwig, MuEnriched =0 or 1, and TeV =7 or 8. Of course this is according to what
   // you've put in Fill_nevent();
   
   m.SetXS();
   // m.SetXS("herwig", 0, 8) ;
   
   m.SetSumXS();
   
   // Here put the name of the PU data root file for the PU reweighting
   TString PUdataFile="MyDataPileupHistogram.root";

   // For the PU reweighting, you can use the 2012 MC distribution with PU sceanrio S7 (SetPU2012_S7(PUdataFile)),
   // with PU scenario S10 (SetPU2012_S10(PUdataFile)), or with your own PU vector. Be carefull this vector has to be
   // of size 60.
   
   m.SetPU2012_S10(PUdataFile);
   
   //The loop on tthe events takes 4 inputs: the trigger, the ptmin of the jets, the ptmax of the jets, and the output
   //file name 
   
   m.Loop(80, 80, 470, "output_data");
   
}
