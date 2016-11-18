TString CampaignLuminosity = "7.65 fb^{-1} (13 TeV)";

TString JECVersionData = "Spring16_25nsV6";
TString JECVersionMC   = "Spring16_25nsV6";

// Eta bins
const int nPtRelEtaBins = 5;
TString PtRelEtaBin[nPtRelEtaBins] = {"anyEta", "Eta06", "Eta12", "Eta18", "Eta24"};
double PtRelEtaEdge[nPtRelEtaBins] = {0., 0.6, 1.2, 1.8, 2.4};

// Pt bins
const int nPtRelPtBins = 18;
TString PtRelPtBin[nPtRelPtBins] = {"Pt2030", "Pt3040", "Pt4050", "Pt5060", "Pt6070", "Pt7080", "Pt80100", "Pt100120", "Pt120140", "Pt140160", "Pt160200", "Pt200260", "Pt260300", "Pt300320", "Pt320400", "Pt400500", "Pt500600", "Pt600"};
double PtRelPtEdge[nPtRelPtBins+1] = {20., 30., 40., 50., 60., 70., 80., 100., 120., 140., 160., 200., 260., 300., 320., 400., 500., 600., 800.};
const int MaxPtRelPtEdge =  800;

double MuonPtCut[nPtRelPtBins] = {5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.};

const int nMaxPU = 50;
const int nMaxPrescales = 4000;

