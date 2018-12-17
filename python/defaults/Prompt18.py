common = {
	'eras' : ['Run2_2018'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/FF9E7432-6F73-A94B-80B4-B2F11D1DCD20.root'],
	'JPCalibration' : 'JPcalib_MC102X_2018_v1',
  'mcGlobalTag' : '102X_upgrade2018_realistic_v12',
	}

data = {
	'inputFiles' : ['/store/data/Run2018A/JetHT/MINIAOD/17Sep2018-v1/00000/00A64001-F644-8740-AC48-14CD4E623E40.root'],	
	'JPCalibration' : 'JPcalib_Data102X_2018_v1',
	'dataGlobalTag' : '102X_dataRun2_Prompt_v11',
}
