// Triggers (and luminosity) info

const int nTriggers = 5;
TString TriggerName[nTriggers]    = {"_DiJet20", "_DiJet40", "_DiJet70", "_DiJet110",   "_Jet300"};
TString JetTriggerName[nTriggers] = {"_PFJet40", "_PFJet40", "_PFJet60",  "_PFJet80", "_PFJet260"};

float MinPtJetTrigger[nTriggers] = { 20.,  50., 100., 140.,  320.};
float MaxPtJetTrigger[nTriggers] = { 50., 100., 140., 320., 1000.};

TString AwayJetTreatment = "LowPtAway";
float PtAwayJet[nTriggers] = {20., 30., 30., 30., 30.};
float PtEmulationTrigger[nTriggers] = {30., 50., 80., 140., 0.};

// Luminosities for the triggers
float Luminosity_DiJet20     =   78759789.113;
float Luminosity_DiJet40     =  113881427.487;
float Luminosity_DiJet70     =  193424037.628;
float Luminosity_DiJet110    =  525115570.604;
float Luminosity_Jet300      = 2525879267.658;

double TriggerLuminosity[nTriggers] = { Luminosity_DiJet20, Luminosity_DiJet40, Luminosity_DiJet70, 
					Luminosity_DiJet110, Luminosity_Jet300};

bool AllowedTrigger[nPtRelPtBins][nTriggers] = { { true, false, false, false, false}, // Pt2030    0
						 { true, false, false, false, false}, // Pt3040    1
						 { true, false, false, false, false}, // Pt4050    2
						 {false,  true, false, false, false}, // Pt5060    3
						 {false,  true, false, false, false}, // Pt6070    4
						 {false,  true, false, false, false}, // Pt7080    5
						 {false,  true, false, false, false}, // Pt80100   6
						 {false, false,  true, false, false}, // Pt100120  7
						 {false, false,  true, false, false}, // Pt120140  8
						 {false, false, false,  true, false}, // Pt140160  9
						 {false, false, false,  true, false}, // Pt160200 10
						 {false, false, false,  true, false}, // Pt200260 11
						 {false, false, false,  true, false}, // Pt260300 12
						 {false, false, false,  true, false}, // Pt300320 13
						 {false, false, false, false,  true}, // Pt320400 14
						 {false, false, false, false,  true}, // Pt400500 15
						 {false, false, false, false,  true}, // Pt500600 16
						 {false, false, false, false,  true} }; // Pt600  17

