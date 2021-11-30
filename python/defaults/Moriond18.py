common = {
	'eras' : ['Run2_2017'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIIFall17MiniAOD/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/60000/002E7FEA-16E0-E711-922D-0242AC130002.root'],
	'JPCalibration' : 'JPcalib_MC94X_2017_v1',
	'mcGlobalTag' : '94X_mc2017_realistic_v12',
	}

data = {
	'inputFiles' : ['/store/data/Run2017C/JetHT/MINIAOD/17Nov2017-v1/20000/00791B22-DCD3-E711-9BF9-001E67396E64.root'],
  'JPCalibration' : 'JPcalib_Data94X_2017_v1',
	'dataGlobalTag' : '94X_dataRun2_ReReco_EOY17_v2',
}
