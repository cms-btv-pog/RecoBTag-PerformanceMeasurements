common = {
	'groups' : ['EventInfo','PV','TagVar','JetInfo','JetSV','JetTrack','CSVTagVar','PFElectron','PFMuon','JetDeepCSV','PatMuon','PatElec','JetDeepFlavour'],
	'eras' : ['Phase2'],
	'miniAOD' : True,
	'usePuppi' : False,
	'usePuppiForFatJets' : True,
	'usePuppiForBTagging' : False,
	'mcGlobalTag' : 'auto:phase2_realistic_T15',
	'remakeAllDiscr' : True,
	'maxJetEta' : 4.5,
	'usePrivateJEC' : True,
	'jecDBFileMC' : 'PhaseIIFall17_V5b_MC',
	'inputFiles' : ['/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/NoPU_110X_mcRun4_realistic_v3-v2/20000/F3F1773B-D087-F941-B4EE-40F14D00A794.root',],
    # 'JPCalibration' : 'JPcalib_MC81X_v0',
}
