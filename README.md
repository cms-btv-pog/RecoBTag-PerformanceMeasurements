# RecoBTag-PerformanceMeasurements

## Software setup

```
cmsrel CMSSW_9_0_3
cd CMSSW_9_0_3/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily

git cms-init

git cms-addpkg DataFormats/BTauReco
git cms-addpkg PhysicsTools/PatAlgos
git cms-addpkg RecoBTag/Combined
git cms-addpkg RecoBTag/CTagging
git cms-addpkg RecoBTag/ImpactParameter
git cms-addpkg RecoBTag/SecondaryVertex
git cms-addpkg RecoBTag/SoftLepton
git cms-addpkg RecoBTau/JetTagComputer
git cms-addpkg RecoVertex/AdaptiveVertexFinder
git cms-addpkg TrackingTools/IPTools

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git

git fetch btv-cmssw NewTaggingVariables_from-CMSSW_9_0_3:NewTaggingVariables_from-CMSSW_9_0_3
git checkout NewTaggingVariables_from-CMSSW_9_0_3

git clone -b 9_0_X_NewVariables --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git
RecoBTag/PerformanceMeasurements

scram b -j8
```

The ntuplizer can be run and configured through

```
cmsRun RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
```

