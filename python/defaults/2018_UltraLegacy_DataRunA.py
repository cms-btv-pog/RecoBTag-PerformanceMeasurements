common = {
	'eras' : ['Run2_2018'],
	'miniAOD' : True,
	'usePrivateJEC': True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIISummer19UL18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/260000/00C28834-56C0-2343-B436-AA8521756E9E.root'],
	'JPCalibration' : 'JPcalib_MC106X_UL2018_v1',
	'mcGlobalTag' : '106X_upgrade2018_realistic_v11_L1v1',
	'jecDBFileMC': 'Summer19UL18_V5_MC',
	}

data = {
	'inputFiles' : ['/store/data/Run2018A/BTagMu/MINIAOD/12Nov2019_UL2018-v1/00000/035456C6-D99D-FD4B-AFA8-DF0B544C65C9.root'],	
    'JPCalibration' : 'JPcalib_Data106X_UL2018_v1',
	'dataGlobalTag' : '106X_dataRun2_v28',
	'jecDBFileData': 'Summer19UL18_RunA_V5_DATA',
}
