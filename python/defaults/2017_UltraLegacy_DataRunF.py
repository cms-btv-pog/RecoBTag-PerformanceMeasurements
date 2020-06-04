common = {
	'eras' : ['Run2_2017'],
	'miniAOD' : True,
	'usePrivateJEC': True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIISummer19UL17MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/40000/FEDB3CA9-7B9C-4D4B-AB90-C4982D57ADEE.root'],
	'JPCalibration' : 'JPcalib_MC106X_UL2017_v1',
	'mcGlobalTag' : '106X_mc2017_realistic_v7',
	'jecDBFileMC': 'Summer19UL17_V5_MC',
	}

data = {
	'inputFiles' : ['/store/data/Run2017D/JetHT/MINIAOD/17Nov2017-v1/20000/0249B143-8CCC-E711-BA7C-0025905C2CD0.root'],	
    'JPCalibration' : 'JPcalib_Data106X_UL2017_v1',
	'dataGlobalTag' : '106X_dataRun2_v28',
	'jecDBFileData': 'Summer19UL17_RunF_V5_DATA',
}
