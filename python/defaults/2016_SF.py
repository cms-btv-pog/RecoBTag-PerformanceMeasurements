common = {
	'eras' : ['Run2_2017'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIISummer16MiniAODv3/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/110000/280DF6F8-7AD0-E811-884A-0025905AC97A.root'],
	'JPCalibration' : 'JPcalib_MC94X_2016_v1',
	'mcGlobalTag' : '94X_mcRun2_asymptotic_v3',
	}

data = {
	'inputFiles' : ['/store/data/Run2016B/JetHT/MINIAOD/17Jul2018_ver2-v2/100000/00323A27-15B8-E811-88B5-E0071B6C9DB0.root'],	
        'JPCalibration' : 'JPcalib_Data94X_2016_v1',
	'dataGlobalTag' : '94X_dataRun2_v10',
}
