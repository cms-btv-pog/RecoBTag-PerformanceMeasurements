from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'BTagAnalyzer_SingleMuonB'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/runBTagAnalyzer_cfg.py'

config.Data.lumiMask = './Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.runRange = '294927-306462'

#SingleMuon
config.Data.inputDataset = '/SingleMuon/Run2017B-17Nov2017-v1/MINIAOD'
#cconfig.Data.inputDataset = '/SingleMuon/Run2017C-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2017D-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2017E-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2017F-17Nov2017-v1/MINIAOD'

#single-top
#config.Data.inputDataset = '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'
#config.Data.inputDataset = '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'

#ttbar
#config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'

#W
#config.Data.inputDataset = '/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v3/MINIAODSIM'
#config.Data.inputDataset = '/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'

#WW
#config.Data.inputDataset = '/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Site.storageSite = 'T3_US_FNALLPC'

