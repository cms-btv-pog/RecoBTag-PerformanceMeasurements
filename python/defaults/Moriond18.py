common = {
	'eras' : ['Run2_2017'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : True,
    'usePrivateJEC': True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIIFall17MiniAOD/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/C8E934F8-1C06-E811-888D-0242AC130002.root'],
	'JPCalibration' : 'JPcalib_MC94X_2017_v1',
	'mcGlobalTag' : '94X_mc2017_realistic_v12',
    'jecDBFileMC': 'Fall17_17Nov2017_V32_94X_MC',
	}

data = {
	'inputFiles' : ['/store/data/Run2017C/JetHT/MINIAOD/17Nov2017-v1/20000/00791B22-DCD3-E711-9BF9-001E67396E64.root'],
  'JPCalibration' : 'JPcalib_Data94X_2017_v1',
	'dataGlobalTag' : '94X_dataRun2_ReReco_EOY17_v2',
    'jecDBFileData': 'Fall17_17Nov2017_V32_94X_DATA',
}
