common = {
	# 'groups' : ['EventInfo','PV','TagVar','JetInfo','JetSV','JetTrack','CSVTagVar','PFElectron','PFMuon','JetDeepCSV','PatMuon','PatElec','JetDeepFlavour'],
	'groups' : ['EventInfo','PV',
                'PFMuon',"PFElectron",
                'PatMuon','PatElec',
                'TagVar','JetTrack',
                'JetInfo','JetSV','CSVTagVar','JetDeepCSV','JetDeepFlavour','CSVTagTrackVar', 'DeepFlavourFeat',],
	'eras' : ['Run3'],
	'miniAOD' : True,
	'usePuppi' : False,
	'usePuppiForFatJets' : True,
	'usePuppiForBTagging' : False,
	# 'dataGlobalTag' : '113X_dataRun3_HLT_v3', #for data
	'mcGlobalTag' : '120X_mcRun3_2021_realistic_v4', #for MC
	'remakeAllDiscr' : True,
	'maxJetEta' : 2.5,
	'usePrivateJEC' : False,
	# 'jecDBFileMC' : 'PhaseIIFall17_V5b_MC',
	# 'inputFiles' : ['/store/relval/CMSSW_11_2_0_pre3/RelValTTbar_14TeV/MINIAODSIM/PU_112X_mcRun3_2021_realistic_v5-v1/20000/4E6AF0AE-EB8E-8847-9C29-5F754974E924.root',],
	'inputFiles' : ['/store/mc/Run3Winter21DRMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/FlatPU30to80_112X_mcRun3_2021_realistic_v16-v2/110000/4edf2114-0dc8-4277-95c5-e55989d35c9e.root',],
    # 'JPCalibration' : 'JPcalib_MC81X_v0',
}
