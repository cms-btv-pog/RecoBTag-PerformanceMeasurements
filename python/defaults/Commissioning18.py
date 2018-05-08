common = {
	'eras' : ['Run2_2018'],
	'miniAOD' : True,
	'storeDeepFlavourTagVariables' : True,
}

mc = {
        #'inputFiles' : ['/store/relval/CMSSW_10_1_0_pre2/RelValTTbar_13/MINIAODSIM/100X_upgrade2018_realistic_v11-v1/20000/14AC524C-981E-E811-B6D1-0CC47A4D7632.root'],
        'inputFiles' : ['/store/mc/RunIISpring18MiniAOD/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v1/30000/008B1124-0B31-E811-BD14-0025904C7A5A.root'],
	'JPCalibration' : 'JPcalib_MC94X_2017_v1',
	'mcGlobalTag' : '94X_mc2017_realistic_v10',
	}

data = {
	'inputFiles' : ['/store/data/Run2018A/JetHT/MINIAOD/PromptReco-v1/000/315/252/00000/04236AA8-404B-E811-9BD8-FA163EF52DAC.root'],
	'JPCalibration' : 'JPcalib_Data94X_2017_v1',
	'dataGlobalTag' : '101X_dataRun2_Prompt_v9',
}
