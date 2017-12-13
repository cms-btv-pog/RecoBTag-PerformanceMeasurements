common = {
	'eras' : ['Run2_2017'],
	'miniAOD' : True,
}

mc = {
	'inputFiles' : ['/store/mc/RunIIFall17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/40000/007E6772-56D6-E711-82B5-0025905A60F8.root'],
	'JPCalibration' : 'JPcalib_MC94X_2017_v1',
	'mcGlobalTag' : '94X_mc2017_realistic_v10',
	}

data = {
	'inputFiles' : ['/store/data/Run2017D/JetHT/MINIAOD/17Nov2017-v1/20000/0249B143-8CCC-E711-BA7C-0025905C2CD0.root'],	
	'JPCalibration' : 'JPcalib_Data94X_2017_v1',
	'dataGlobalTag' : '94X_dataRun2_ReReco_EOY17_v2',
}
