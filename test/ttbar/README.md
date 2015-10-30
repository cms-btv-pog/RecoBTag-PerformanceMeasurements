----------------------------------------------
# RecoBTag-PerformanceMeasurements with ttbar

### Installation
See base installation here
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagAnalyzer#Recipe_for_7_4_X_releases

### Producing the trees
```
cmsRun runBTagAnalyzer_cfg.py runOnData=False miniAOD=True useTTbarFilter=True maxEvents=-1
```
Will run locally the analyzer for testing purposes
```
python submitToGrid.py -j data/samples_Run2015_25ns.json -c ${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py -s
```
Will submit the samples described in the json to the grid.
Partial submission can be made with -o csv_list
Don't forget to init the environment for crab3
(e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial)
Other json files with samples (e.g. for systematics scans) are also available under data.
```
python checkProductionIntegrity.py -i /store/group/phys_btag/performance/TTbar/7b810a5
```
Can run as soon as ntuple production starts to end, to move from crab output directories to a more simple directory structure
which can be easily parsed by the local analysis. 
Once production is finished you can call with --cleanup to remove original crab directories in EOS.

### Running local analysis
```
python runTTbarAnalysis.py -i /store/group/phys_btag/performance/TTbar/7b810a5 -j ttbar_Run2015.json -n 8
```
Once grid jobs are run, and ntuples are stored in a given directory, you can run the local analysis to produce the slimmed ntuples for the efficiency measurement.
MC will be weighted by cross section. The number after -n indicates how many threads should be used.
NB The first time the script will produce a pickle file with the weights to be used according to the number of files found, xsec specified in the json file.
Yields are normalized to 1/pb of integrated luminosity.
The integrated luminosity has to be computed using the json reported by crab using the directories under grid.
See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial and https://twiki.cern.ch/twiki/bin/view/CMS/Lcr2 on how to do these intermediate operations
It prints out s.th. like "Produced normalization cache (analysis/.xsecweights.pck)"
In case you update the trees, xsec or lumi you have to remove by hand the pickle file.
```
python plotter.py -i analysis/ -j ttbar_Run2015.json  -l 71.44
python plotter.py -i analysis/ -j ttbar_Run2015.json  -l 71.44 --only evsel --saveTeX
```
Makes control plots and stores all in a ROOT file. Different options may be passed to filter plots, and show differently the plots. 
```
python fitDYFromData.py -i analysis/plots/plotter.root
```
This script uses the output plots for MET to fit a common scale factor for DY. MET is modelled from under the Z peak
and used to fit the MET observed after the requirement of at least two jets.
The scale factor can be applied to all channels as it reflects possible mismodelling in the jet multiplicity of the DY MC.

### Performance analyses

#### Templated methods
```
sh KIN_runClassifier.sh
```
After running the local analysis use the kin tree stored in the ttbar sample to train a kinematics discriminator for b-jets in ttbar events.
The script compiles and runs KIN_trainClassifier.C which should be modified in case different trainings are required.
```
python runTTbarAnalysis.py -i /store/group/phys_btag/performance/TTbar/7b810a5 -j ttbar_Run2015.json --tmvaWgts analysis/KIN_weights/ --dyScale analysis/plots/.dyScaleFactor.pck  -n 8
```
Re-run the analysis to store the KIN discriminator value per jet
```
python Templated_btagEffFitter.py -i analysis/ -o analysis_inc/ -t taggers_Run2015.json -n 8 
python Templated_btagEffFitter.py -i analysis/ -o analysis_ll/ -t taggers_Run2015.json -n 8 --channels -121,-169
python Templated_btagEffFitter.py -i analysis/ -o analysis_emu/ -t taggers_Run2015.json -n 8 --channels -143
```
Runs the fits to the templates to determine the scale factors. Valid for KIN, Mlj, JP, others one may wish to add.
The base procedure is similar for all. The first time to run will take a long time as templates need to be created.
If templates are stable and only fits need to be redone when can run with the option "--recycleTemplates"
```
python createSFbSummaryReport.py -i analysis/kindisc_templates/.csv_fits.pck -o csv_fits.tex
```
Parses the fit results and creates a TeX file with the tables as well as the plots with the measured efficiencies.
```
python runUnfoldDiscriminatorShape.py -i analysis_emu/kindisc_templates/csv.root -o analysis_emu/kindisc_templates/
```
Runs unbinned fits to the data and reweights the data using the KIN discriminator to obtain an estimate
of the b-tagging discriminator shape. 


#### FtM method