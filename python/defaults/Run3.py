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
    #'inputFiles' : ['/store/mc/Run3Winter21DRMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/FlatPU30to80_112X_mcRun3_2021_realistic_v16-v2/110000/4edf2114-0dc8-4277-95c5-e55989d35c9e.root',],
    #'inputFiles' :['/store/relval/CMSSW_12_1_0_pre2/RelValTTbar_14TeV/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v3/10000/12e01fdc-f329-4d40-8c75-c5dde3ee6e5f.root'],
    #'inputFiles' :['/store/relval/CMSSW_12_1_0_pre2/RelValTTbar_14TeV/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v1_RECOonly-v1/10000/78657524-f827-4856-8e43-30264f354740.root'],
    #'inputFiles' :['/store/relval/CMSSW_12_1_0_pre2/RelValTTbar_14TeV/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v1-v1/10000/73d68dcd-591d-45fc-a58b-eeddd7b886d1.root'],

    #'inputFiles' : ['/store/relval/CMSSW_12_1_0_pre2/RelValTTbar_14TeV/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/33c870c2-b7d8-437f-aea1-74b1e0629783.root'],
    #'outFilename' : "JetTree_mkFit_NoPU.root",

    'inputFiles' : ['/store/relval/CMSSW_12_1_0_pre2/RelValTTbar_14TeV/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/78d04134-4010-4c4d-9ed2-6cf597f0f554.root'],
    'outFilename' : "JetTree_ref_NoPU.root",

    # 'JPCalibration' : 'JPcalib_MC81X_v0',
}
