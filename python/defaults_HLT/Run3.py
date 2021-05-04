common = {
	'groups' : ['EventInfo','PV','PFMuon',"PFElectron",
                'TagVar','JetInfo','JetSV','JetTrack','CSVTagVar','JetDeepCSV','JetDeepFlavour','CSVTagTrackVar', 'DeepFlavourFeat'
                'PuppiJetInfo','PuppiJetSV','PuppiJetTrack','PuppiJetCSVTagVar','PuppiJetDeepCSV','PuppiJetDeepFlavour','PuppiJetTagVar','PuppiJetCSVTagTrackVar','PuppiJetDeepFlavourFeat',
                'CaloJetInfo','CaloJetSV','CaloJetTrack','CaloJetCSVTagVar','CaloJetDeepCSV','CaloJetTagVar','CaloJetCSVTagTrackVar',
                ],
	'eras' : ['Run3'],
    'runCaloJetVariables' : True,
    'runPuppiJetVariables' : True,
	'mcGlobalTag' : '112X_mcRun3_2021_realistic_v14',
	'maxJetEta' : 4.5,
	'usePrivateJEC' : False,
	'jecDBFileMC' : 'PhaseIIFall17_V5b_MC',
	'inputFiles' : ['/store/relval/CMSSW_11_2_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v5-v1/20000/0F40986A-FECA-F04A-99E9-4F35A32C369E.root',],
    # 'JPCalibration' : 'JPcalib_MC81X_v0',
}
