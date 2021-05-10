common = {
	'eras' : ['Run2_2016'],
	'miniAOD' : True,
	'usePrivateJEC': True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIISummer20UL16MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v2/00000/011B6125-7815-DA41-8156-3AE937358B20.root'],
	'JPCalibration' : 'JPcalib_MC106X_UL2016_postVFP_v1',
	'mcGlobalTag' : '106X_mcRun2_asymptotic_v15',
	'jecDBFileMC': 'Summer19UL16_V7_MC',
	}

data = {
	'inputFiles' : ['/store/data/Run2016G/BTagMu/MINIAOD/21Feb2020_UL2016-v1/20000/085E85E8-4F7B-994E-8119-4591FC4D44B7.root'],	
    'JPCalibration' : 'JPcalib_Data106X_UL2016_v1',
	'dataGlobalTag' : '106X_dataRun2_v32',
	'jecDBFileData' : 'Summer19UL16_RunBCDEFGH_Combined_V7_DATA',
}
