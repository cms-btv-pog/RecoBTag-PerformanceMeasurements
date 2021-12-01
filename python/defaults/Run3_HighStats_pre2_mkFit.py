common = {
	# 'groups' : ['EventInfo','PV','TagVar','JetInfo','JetSV','JetTrack','CSVTagVar','PFElectron','PFMuon','JetDeepCSV','PatMuon','PatElec','JetDeepFlavour'],
	'groups' : ['EventInfo','PV',
                'PFMuon',#"PFElectron",
                'PatMuon','PatElec',
                'TagVar','JetTrack',
                'JetInfo','JetSV','CSVTagVar','JetDeepCSV','JetDeepFlavour','CSVTagTrackVar', 'DeepFlavourFeat',],
	'eras' : ['Run3'],
	'miniAOD' : True,
	'usePuppi' : False,
	'usePuppiForFatJets' : True,
	'usePuppiForBTagging' : False,
	# 'dataGlobalTag' : '113X_dataRun3_HLT_v3', #for data
    #'mcGlobalTag' : '120X_mcRun3_2021_realistic_v4', #for MC
    'mcGlobalTag': '122X_mcRun3_2021_realistic_v1',
	'remakeAllDiscr' : False,
	'maxJetEta' : 2.5,
	'usePrivateJEC' : False,
	# 'jecDBFileMC' : 'PhaseIIFall17_V5b_MC',

    #dasgoclient --query="file dataset=/RelValTTbar_14TeV/CMSSW_12_2_0_pre2-122X_mcRun3_2021_realistic_v1_TkmkFitHighStat-v1/MINIAODSIM"
    'inputFiles' : [
        '/store/relval/CMSSW_12_2_0_pre2/RelValTTbar_14TeV/MINIAODSIM/122X_mcRun3_2021_realistic_v1_TkmkFitHighStat-v1/2580000/2ce11d10-79ae-4ec7-9c78-93ee4a84f112.root'
    ],
    'outFilename' : "JetTree_highStats_pre2_mkFit.root",

    # 'JPCalibration' : 'JPcalib_MC81X_v0',
}
