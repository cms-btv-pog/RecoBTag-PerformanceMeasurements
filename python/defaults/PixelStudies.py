common = {
	'eras' : ['Run2_2017'],
	'miniAOD' : True,
	'runDeepFlavourTagVariables' : True,
}

mc = {
	'inputFiles' : ['/store/relval/CMSSW_10_5_0_pre1/RelValTTbar_13/MINIAODSIM/104X_upgrade2018_realistic_pixelCand_v1-v1/20000/BE6E5F4F-4AAE-C94F-BB04-946DC9406DB8.root'],
	'JPCalibration' : 'JPcalib_MC94X_2017_v1',
	'mcGlobalTag' : '104X_upgrade2018_realistic_pixelCand_v1',
	}

data = {
	'inputFiles' : ['/store/relval/CMSSW_10_5_0_pre1/JetHT/MINIAOD/104X_dataRun2_pixel2DTemplate_Cand_v1_ClRepairOn_RelVal_jetHT2018D-v1/20000/F8117214-B35A-BF42-8FEF-48CB817731C0.root'],
  'JPCalibration' : 'JPcalib_Data94X_2017_v1',
	'dataGlobalTag' : '104X_dataRun2_pixel2DTemplate_Cand_v1',
}