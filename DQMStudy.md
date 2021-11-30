## Following are instructions for detailed DQM studies :

Setup the RecoBTag-PerformanceMeasurements pactkage as described in the README.md

# Make Target Reference BTagAnalyzer NTuple
```
cp RecoBTag/PerformanceMeasurements/python/defaults/Run3.py RecoBTag/PerformanceMeasurements/python/defaults/Run3Reference.py 
[ edit Run3Refernce to point to the right file ]
cd RecoBTag/PerformanceMeasurements/test/
cmsRun runBTagAnalyzer_cfg.py defaults=Run3Reference runOnData=(True or False, depending on your needs) 
cd -
```

# Make Target Reference BTagAnalyzer NTuple
```
cp RecoBTag/PerformanceMeasurements/python/defaults/Run3.py RecoBTag/PerformanceMeasurements/python/defaults/Run3Target.py 
[ edit Run3Target to point to the right file ]
cd RecoBTag/PerformanceMeasurements/test/
cmsRun runBTagAnalyzer_cfg.py defaults=Run3Target runOnData=(True or False, depending on your needs) 
cd -
```

### Process Ntuples To make ROOT file with histograms

Additional Software setup to process BTagNtuples for offline validation monitoring

```
git cms-merge-topic johnalison:MakePy3Bind11ParameterSetsIncludingCommandLineArguments
git clone -b BTVDQM git@github.com:patrickbryant/nTupleAnalysis.git
git clone -b 12_1_X git@github.com:johnalison/TriggerStudies.git
```

```
BTagAnalyzer TriggerStudies/NtupleAna/scripts/BTagAnalyzer_cfg.py \
    --inputRAW  path/to/target/BTagAnalyzer/Ntuple/JetTree.root \
    --inputAOD path/to/reference/BTagAnalyzer/Ntuple/JetTree.root \
    -o $PWD \
    -y Run3 \
    --histogramming 1 \
    --histFile hists_Run3_offlineValidation.root \
    --isMC \
    --doTracks \
    --jetDetailString matchedJet.allTaggers.innerPixHit.noInnerPixHit.DeepJetBins \
    --pfJetName "" \ 
    --nevents -1
```



### Make plots presentation
``` 
git clone -b 12_1_X  git@github.com:johnalison/ROOTHelp.git
source ROOTHelp/setup.sh
source TriggerStudies/plotting/offlineDQMPlots.sh hists_Run3_offlineValidation.root hists_Run3_offlineValidation targetName [python3] 
```
Outpus a slides in hists_Run3_offlineValidation directory
give the option python3 to actually make the plots 
