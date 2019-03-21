common = {
	'eras' : ['Run2_2018'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : True,
        'doBoostedCommissioning': True,
        'groups': ['EventInfo', 'Devdatta', 'DoubleBCommissioning'],
}

mc = {
	'inputFiles' : ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext3-v2/80000/26E977A9-0835-9F43-BA5C-A0CE1E452C0A.root'],
	'JPCalibration' : 'JPcalib_MC102X_2018_v1',
        'mcGlobalTag' : '102X_upgrade2018_realistic_v18',
	}

data = {
	'inputFiles' : ['/store/data/Run2018A/JetHT/MINIAOD/17Sep2018-v1/00000/00A64001-F644-8740-AC48-14CD4E623E40.root'],
	'JPCalibration' : 'JPcalib_Data102X_2018_v1',
	'dataGlobalTag' : '102X_dataRun2_Prompt_v13',
}
