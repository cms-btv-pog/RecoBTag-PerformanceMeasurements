#!/bin/csh

foreach line ( "`cat crabConfig.list`" )

   #ignore line if it is commented out
   set x = `echo ${line} | cut -c 1` 
   if ("${x}" == "#") then
      continue
   endif
   
   #set variables
   set requestName_ = `echo ${line} | awk '{print $1}'`
   echo "running on ${requestName_}"
   set inputDataset_ = `echo ${line} | awk '{print $2}'`
   set isData_ = `echo ${line} | awk '{print $3}'`
   set crabConfig_ = "./crab_projects/crabConfig_${requestName_}.py"
   
   # write crabConfig file
   echo "from CRABClient.UserUtilities import config, getUsernameFromSiteDB" > ${crabConfig_}
   echo "config = config()" >> ${crabConfig_}
   echo " " >> ${crabConfig_}
   echo "config.General.requestName = 'BTagAnalyzer_${requestName_}'" >> ${crabConfig_}
   echo "config.General.workArea = 'crab_projects'" >> ${crabConfig_}
   echo "config.General.transferOutputs = True" >> ${crabConfig_}
   echo "config.General.transferLogs = True" >> ${crabConfig_}
   echo " " >> ${crabConfig_}
   echo "config.JobType.pluginName = 'Analysis'" >> ${crabConfig_}
   echo "config.JobType.psetName = '../test/runBTagAnalyzer_cfg.py'" >> ${crabConfig_}
   echo "config.JobType.pyCfgParams = ['runOnData=${isData_}']" >> ${crabConfig_}
   echo " " >> ${crabConfig_}
   if ($isData_ == "True") then
      echo "config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'" >> ${crabConfig_}
      echo "config.Data.runRange = '294927-306462'" >> ${crabConfig_}
   endif
   echo "config.Data.inputDataset = '${inputDataset_}'" >> ${crabConfig_}
   echo "config.Data.splitting = 'FileBased'" >> ${crabConfig_}
   echo "config.Data.unitsPerJob = 1" >> ${crabConfig_}
   echo " " >> ${crabConfig_}
   echo "config.Site.storageSite = 'T3_US_FNALLPC'" >> ${crabConfig_}

   # submit crabConfig file
   crab submit --config ${crabConfig_}
end
