common = {
	'eras' : ['Run2_2016'],
	'miniAOD' : True,
	'usePrivateJEC': True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIISummer20UL16MiniAODAPV/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v8-v2/00000/00098436-82CB-8747-AB1A-96CDB0A9B640.root'],
	'JPCalibration' : 'JPcalib_MC106X_UL2016_v1',
	'mcGlobalTag' : '106X_mcRun2_asymptotic_preVFP_v9',
	'jecDBFileMC': 'Summer19UL16APV_V7_MC.db',
	}

data = {
	'inputFiles' : [ '/store/data/Run2016C/BTagMu/MINIAOD/21Feb2020_UL2016_HIPM-v1/10000/1BA24193-AD07-7849-8595-0EC1FEC99450.root'],	
    'JPCalibration' : 'JPcalib_Data106X_UL2016_v1',
	'dataGlobalTag' : '106X_dataRun2_v32',
	'jecDBFileData' : 'Summer19UL16_RunBCDEFGH_Combined_V7_DATA.db',
}
