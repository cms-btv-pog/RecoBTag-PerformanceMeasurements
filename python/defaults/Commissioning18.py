common = {
	'eras' : ['Run2_2018'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : True,
}

mc = {
        'inputFiles' : ['/store/mc/RunIISpring18MiniAOD/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v1/910000/80881E0F-4838-E811-9E65-485B39897242.root '],
	'JPCalibration' : 'JPcalib_MC94X_2017_v1',
	'mcGlobalTag' : '94X_mc2017_realistic_v12',
	}

data = {
	'inputFiles' : ['/store/data/Run2018A/JetHT/MINIAOD/PromptReco-v1/000/315/252/00000/04236AA8-404B-E811-9BD8-FA163EF52DAC.root'],
        'JPCalibration' : 'JPcalib_Data94X_2017_v1',
	'dataGlobalTag' : '101X_dataRun2_Prompt_v11',
}
