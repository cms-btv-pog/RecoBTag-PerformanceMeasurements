# To calculate the bb scale factors for W jets in ttbar events  

## 1) create ntuples  
a) Add a simple EDFilter to drastically reduce the event load:  
https://github.com/fojensen/RecoBTag-PerformanceMeasurements/blob/9_4_X/plugins/recoFilter.cc  
https://github.com/fojensen/RecoBTag-PerformanceMeasurements/blob/9_4_X/test/runBTagAnalyzer_cfg.py#L1541  
b) add an EDAnalyzer to fill getTrueNumInteractions for mc pileup corrections:  
https://github.com/fojensen/RecoBTag-PerformanceMeasurements/blob/9_4_X/plugins/plotPUtrue.cc  
https://github.com/fojensen/RecoBTag-PerformanceMeasurements/blob/9_4_X/test/runBTagAnalyzer_cfg.py#L1540  
c) submit jobs to crab: `./crabConfig.sh`  

## 2) make pu weights  
a) `root pu.c+`  

## 3) make skims  
a) `root makeskims.c+` (_each dataset needs its own interactive job, in the future a script should be written to submit the jobs to Condor in parallel_)  
b) combine SingleMuon skims: `hadd ./skims/SingleMuon.root ./skims/SingleMuon*.root`  

## 4) make data/mc comparisons  
a) `root datamccomparisons.c+`  

## 5) make scale factors  
a) `root sf.c+`  
b) `root sf_systematics.c+`  
c) combine the different wp sf root files to one file and add hard-coded systematic errors: `root sfaddfiles.root+`  
