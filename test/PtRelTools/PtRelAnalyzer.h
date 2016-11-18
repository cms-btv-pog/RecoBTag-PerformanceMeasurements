#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"

#include "CampaignParameters/Run201680X-ICHEP2016/BaseParameters.h"
#include "CampaignParameters/Run201680X-ICHEP2016/TriggerInfo.h"
#include "CampaignParameters/Run201680X-ICHEP2016/Taggers.h"
#include "CampaignParameters/Run201680X-ICHEP2016/Systematics.h"

// Choose the production 
#include "CampaignParameters/Run201680X-ICHEP2016/BaseProduction.h"
//#include "CampaignParameters/Run201680X/Run2015BProduction.h"
  
float TotalScaleFactorSystematic[nFitPtBins][nPtRelEtaBins];
float ScaleFactorSystematic[nFitPtBins][nPtRelEtaBins][nScaleFactorSystematics];
float ScaleFactorValue[nFitPtBins][nPtRelEtaBins];

struct ScaleFactorResult {double EffData; double EffDataRaw; double EffDataError; double EffMC; double EffMCError; double SF; double SFError; double FracTag; double FracTagError; double FracUntag; double FracUntagError; double Chi2NTag; double Chi2NUntag; int DataEventTag; int DataEventUntag; };

struct PtRelFitResult {double bFraction; double bFractionError; double cFraction; double cFractionError; double lightFraction; double lightFractionError; double Status; double Chi2N;};

const int nAwayTaggers = 17;
TString AwayTaggerName[nAwayTaggers] = {"TCHPM", "TCHPL", "TCHPT", "JBPT", "JBPT1", "JBPM", "JBPM1", "JBPM2", "JBPM3", "JBPM4", "JBPM5", "JBPM6", "JBPM7", "JBPM8", "JBPM9", "JBPL", "NONE"};

class PtRelAnalyzer {
  
 public:
  
  PtRelAnalyzer(TString TemplateFlag, TString PileUpReweighting = "", TString KinematicWeighting = "", TString SelectionFlag = "");
  ~PtRelAnalyzer();

  void BookHistograms();
  void FillHistograms(TString DataType, int DataRange);
  void FillSubjetHistograms(TString DataType, int DataRange);
  void MergeHistograms(TString DataType, int FirstBin = 0, int LastBin = 100);

  void BookLightHistograms();
  void FillLightHistograms(TString DataType, int DataRange);
  void FillSubjetLightHistograms(TString DataType, int DataRange);
  void MergeLightHistograms(TString DataType, int FirstBin = 0, int LastBin = 100);

  void BuildTemplates(bool AddLightTemplates = false, TString DataFlag = "BTagMu", TString MCFlag = "QCDMu");
  void BuildTemplatesEtaSplit(bool AddLightTemplates = false);

  void BookSystem8Histograms();
  void FillSystem8Histograms(TString DataType, int DataRange);
  void MergeSystem8Histograms(TString DataType, int FirstBin = 0, int LastBin = 100);

  void BuildSystem8Templates();
  
  void ComputeKinematicWeights(TString SystematicFlag = "_Central", TString DataType = "QCDMu", TString LightTemplates = "");

  void ComputePtRelScaleFactors(TString Taggers, TString SystematicFlag = "_Central", TString FitOption = "_LightTemplatesRatio", TString PrintPlot = "png", TString PtBin = "All", TString EtaBin = "anyEta", TString TemplateFlag = "All");

  void PlotBTagPerformance(TString PlotFlag, TString Tagger, TString EtaList[], TString ConfigurationList[], TString SystematicList[], float MaxJetPt = 1000000., TString DrawedSystematics = "NONE", TString Type = "Performance");

  void ComputeScaleFactorSystematics(TString Taggers = "All", TString EtaBin = "anyEta", TString CentralFlag = "_Central");

  void EventCounter(TString DataType, int DataRange);
  
  void ComputePileUpWeights(TString PUMethod, TString DataType, int DataRange = -1, TString PrimaryDataset = "BTagMu");
  
  void CompareDataToMC(TString Variable, TString EtaBin = "anyEta", TString PtBin = "All", int Rebin = 1, TString LightTemplates = "All", TString DataType = "BTagMu");
 
  void CompareTemplates(TString Variable, TString FirstFileName, TString SecondFileName, TString FirstName, TString SecondName);

  void StudyAwayJetSelection(TString Step, int DataRange = -1);

  void ComputeMethodsOverlap();

  void ComputeBTaggingEfficiency(TString DataType, int DataRange = -1);
 
 private:

  TString TemplateVariable, PUWeighting, KinWeighting, Selection;

  int nBinsForTemp, TemplateRebinning;
  float LowerEdgeForTemp, UpperEdgeForTemp;
  TString TemplateXTitle, TemplateYTitle;

  int GetAwayJet(TString ThisAwayTaggerName, int jMu, float AwayDRCut, bool IsUnique);
  int GetAwayStdJet(TString ThisAwayTaggerName, int jMu, float AwayDRCut, bool IsUnique, int TriggerIdx);
  void GetDRCuts(TString ThisSystematicName, float JetPt, float *TrackDRCut, float *TrackMinDRCut);
  
  double KinematicWeight[nTriggers][nSystematics][MaxPtRelPtEdge][60][48];
  float PileUpWeight[nTriggers][nMaxPU][nSystematics];
  float BTemplateCorrections[100][nPtRelPtBins][2];
  int TriggerPrescaleRunNumber[nTriggers][nMaxPrescales];
  int TriggerPrescaleLumiSection[nTriggers][nMaxPrescales];
  int TriggerPrescale[nTriggers][nMaxPrescales];

  void ResetHistograms(TString DataType);
  void SaveHistograms(TString OutputFileName, TString DataType);
  void NormalizeHistogramsToEntries();

  void ResetLightHistograms(TString DataType);
  void SaveLightHistograms(TString OutputFileName);
  void NormalizeLightHistogramsToEntries();

  void ResetSystem8Histograms(TString DataType);
  void SaveSystem8Histograms(TString OutputFileName, TString DataType);
  void NormalizeSystem8HistogramsToEntries();

  void GetKinematicWeights(TString DataType);
  void GetPileUpWeights(TString DataType);
  void GetBTemplateCorrections();
  void GetTriggerPrescales(TString DataType);
  float TriggerPrescaleWeight(TString DataType, int TriggerIndex);

  TString HistogramName(TString VariableName, TString PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx);
  TString HistogramName(TString VariableName, int PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx);

  void BookTemplates();
  void SaveTemplates(TString TemplateFileName, bool AddLightTemplates);

  void BookSystem8Templates();
  void SaveSystem8Templates(TString TemplateFileName);

  bool PassTriggerEmulation(int TriggerIdx, int MuonJetIdx);
  bool PassEventTriggerEmulation(int TriggerIdx, TString DataType);
  
  // Histograms
  TH1D *MuDRForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *MuPtForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *JetPtForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *JetEtaForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PVMultiplicity[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  
  TH1D *PtRelTagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PtRelUntagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PtRelLightTagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *PtRelLightUntagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
 
  TH1D *PtRelLightTagForSystematic[nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PtRelLightUntagForSystematic[nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  
  TH2D *Observed[nPtRelPtBins][nPtRelEtaBins][nSystematics];

  // Templates
  TH1D *JetPt[4][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *JetEta[4][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *MuonPt[2][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *MuonDR[2][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PVEvent[2][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *PtRel[2*nTaggers][nSystematics][nPtRelPtBins][nPtRelEtaBins][8]; 

  TH1D *System8[2*nTaggers+2][nSystematics][nPtRelPtBins][nPtRelEtaBins][3]; 

  TH1D *System8ForDataWeighting[2*nTaggers+2][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *System8ForBJetWeighting[2*nTaggers+2][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *System8ForLJetWeighting[2*nTaggers+2][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];

  void WriteKinematicWeights(TH1D *&HistoDivide, TString Variable, int PtBin, int EtaBin, TString DataType);

  TFile *PtRelTemplateFile;
  TH1D *GetPtRelTemplate(bool IsMC, TString Flavour, TString Tagger, TString EtaBin, TString PtBin, TString FitOption);
  PtRelFitResult PtRelFit(TH1D* DataTemplate, TH1D* bTemplate, TH1D* cTemplate, TH1D* lgTemplate, TString FitOption);
  double BinomialError(double a, double b);
  double EfficiencyStatisticalError(double a, double b);
  double PtRelEfficiencyError(double a, double b, double fa, double fb, double e_fa, double e_fb);
  double RatioPropagationError(double a, double b, double e_a, double e_b);
  ScaleFactorResult PtRelScaleFactor(TString Tagger, TString EtaBin, TString PtBin, TString SystematicFlag, TString FitOption, TString PrintPlot);

  ScaleFactorResult ScaleFactorResultFromTable(TString TableName);
  TString ScaleFactorResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString KinWeightingFlag, TString SelectionFlag, TString ProductionFlag);
  TString ScaleFactorResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString Configuration, TString ProductionFlag);

  TCanvas *BTagPerformanceCanvas(TString Type, float CanvasHeight = 400., TString LogAxes = "");
  TH1D *MCEfficiency, *DataEfficiency, *DataMCSF, *DataMCSFSystematics;
  void FillBTagPerformanceHistograms(TString Tagger, TString EtaBin, TString Configuration, TString Systematic, int ColorIdx, TString DrawedSystematics);

  void ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString CentralFlag);
  void ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString Configuration, TString CentralFlag);
  void ComputeScaleFactorSystematics(TString Tagger, TString EtaBin, TString Configuration, TString CentralFlag);

  void CompareHistograms(TCanvas *PlotCanvas, TString Variable, TH1D *FirstHisto, TH1D *SecondHisto, TString PlotName, int Rebin = 1);

};

