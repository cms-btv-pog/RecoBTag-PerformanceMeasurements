common = {
	'eras' : ['Run2_2016'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : False,
    'doBoostedCommissioning': True,
    'groups': ['EventInfo', 'Devdatta', 'DoubleBCommissioning'],
	'usePrivateJEC': True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIISummer16MiniAODv3/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/120000/92DEC622-A2E6-E811-9DE8-0025904CDDF8.root',],
	'JPCalibration' : 'JPcalib_MC94X_2016_v1',
	'mcGlobalTag' : '94X_mcRun2_asymptotic_v3',
	'jecDBFileMC': 'Summer16_07Aug2017_V11_MC',  ### with the FormulaEvaluator bug
	}

data = {
	'inputFiles' : ['/store/data/Run2016B/JetHT/MINIAOD/17Jul2018_ver2-v2/100000/00323A27-15B8-E811-88B5-E0071B6C9DB0.root'],
        'JPCalibration' : 'JPcalib_Data94X_2016_v1',
	'dataGlobalTag' : '94X_dataRun2_v10',
	'jecDBFileData': 'Summer16_07Aug2017All_V11_DATA',
}
