
bool VerboseSystematics = false;

//const int nSystematics = 1; 
//TString SystematicName[nSystematics] = {"_Central"};

//const int nSystematics = 2; 
//TString SystematicName[nSystematics] = {"_Central", "_JBPL"};

//const int nSystematics = 15; 
//TString SystematicName[nSystematics] = {"_Central", "_GluonSplittingBDown", "_GluonSplittingBUp", "_MuPt6", "_MuPt8", "_MuDRMinus", "_MuDRPlus", "_TCHPL", "_TCHPT", "_JBPL", "_JBPT", "_TCHPM", "_JBPM7", "_BTemplatesMinus", "_BTemplatesPlus"};

const int nSystematics = 15;
TString SystematicName[nSystematics] = {"_Central", "_PileUpDown", "_PileUpUp", "_GluonSplittingBDown", "_GluonSplittingBUp", "_MuPt6", "_MuPt8", "_MuDRMinus", "_MuDRPlus", "_JEUDown", "_JEUUp", "_JBPL", "_JBPT", "_BTemplatesMinus", "_BTemplatesPlus"};

const int nFitSystematics = 19;  
//TString FitSystematicName[nFitSystematics] = {"_Central", "_GluonSplittingBDown", "_GluonSplittingBUp", "_MuPt6", "_MuPt8", "_MuDRMinus", "_MuDRPlus", "_JEUDown", "_JEUUp", "_BTemplatesMinus", "_BTemplatesPlus", "_Central_LightCharmRatio-m30", "_Central_LightCharmRatio-p30", "_Central_bTempRatioCorr", "_Central_bTempRatioCorr", "_Central_lCorr", "_Central_lCorr"}; 
TString FitSystematicName[nFitSystematics] = {"_Central", "_PileUpDown", "_PileUpUp", "_GluonSplittingBDown", "_GluonSplittingBUp", "_MuPt6", "_MuPt8", "_MuDRMinus", "_MuDRPlus", "_JEUDataDown", "_JEUDataUp", "_JBPL", "_JBPT", "_BTemplatesMinus", "_BTemplatesPlus", "_Central_LightCharmRatio-m30", "_Central_LightCharmRatio-p30", "_Central_bCorr", "_Central_bCorr"}; 

const int nScaleFactorSystematics = 10;
TString ScaleFactorSystematicName[nScaleFactorSystematics] = {"_Statistics", "_PileUp", "_GluonSplitting", "_MuPt", "_MuDR", "_JEUData", "_AwayTagger", "_BFragmentation", "_LightToCharmRatio", "_bCorr"};
TString ScaleFactorSystematicStoredName[nScaleFactorSystematics] = {"statistic", "pileup", "gluonsplitting", "mupt", "mudr", "jes", "jetaway", "bfragmentation", "l2c", "btempcorr"};


