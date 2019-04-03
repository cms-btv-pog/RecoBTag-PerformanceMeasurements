common = {
	'eras' : ['Run2_2018'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : False,
        'doBoostedCommissioning': True,
        'groups': ['EventInfo', 'Devdatta', 'DoubleBCommissioning'],
}

mc = {
	'inputFiles' : ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/PUAvg50ForMUOVal_102X_upgrade2018_realistic_v15_ext1-v1/80000/1410E06E-D901-3448-9449-AB10C942CE03.root'],
	'JPCalibration' : 'JPcalib_MC102X_2018_v1',
	'mcGlobalTag' : '102X_upgrade2018_realistic_v12',
  'usePrivateJEC': True,
  'jecDBFileMC': 'Autumn18_V8_MC'
	}

####17Sep2018
#data = {
#	'inputFiles' : ['/store/data/Run2018C/BTagMu/MINIAOD/17Sep2018-v1/60000/FE0C866E-6D00-DE47-89AB-A15CC16813B9.root'],
# 'JPCalibration' : 'JPcalib_Data102X_2018_v1',
#	'dataGlobalTag' : '102X_dataRun2_Sep2018Rereco_v1',
# 'usePrivateJEC': True,
# 'jecDBFileData': 'Autumn18_RunABCD_V8_DATA'
#}

###Prompt
data = {
	'inputFiles' : ['/store/data/Run2018D/BTagMu/MINIAOD/PromptReco-v2/000/325/022/00000/6F0A55A7-7F5F-A54D-A01E-243B27BBCB9A.root'],
  'JPCalibration' : 'JPcalib_Data102X_2018_v1',
	'dataGlobalTag' : '102X_dataRun2_Prompt_v11',
}
