#########################
# KIN BTag Method Notes
# Moriond 2018 Campaign
#########################
# Joshuha Thomas-Wilsker
#  CERN, IHEP Beijing
#########################

CMSSW to use for new campaign: CMSSW_9_4_1

*****

Fall17 Production:
#das_client.py --query="dataset dataset=/*/RunIIFall17MiniAOD-94X_mc2017_realistic_v*/MINIAODSIM status=*" --limit=300

*****

Triggers:
https://login.cern.ch/adfs/ls/?wa=wsignin1.0&wreply=https%3A%2F%2Ftwiki.cern.ch%2FShibboleth.sso%2FADFS&wct=2017-12-12T11%3A50%3A10Z&wtrealm=https%3A%2F%2Ftwiki.cern.ch%2FShibboleth.sso%2FADFS&wctx=cookie%3A1513079410_64ae#Year_2017

MC Samples:
Many of the samples have changed in various ways. For example, for ttbar samples we now have to use the samples where the ttbar events are divided up according to the decay of the W bosons. This means when we include the samples in the samples.json we need to calculate the cross-section for the process. The was done for the ttbar samples by taking the ttbar production XS from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO#Top_quark_pair_cross_sections_at
This was then multiplied by the branching fractions of the W boson taken from the PDG: http://pdg.lbl.gov/2017/tables/rpp2017-sum-gauge-higgs-bosons.pdf
As for the Drell Yan samples, they are divided into HT slices which means the sample should have a different cross-section depending on the HT slice (unless we can use the FXFX sample and then it’s inclusive across HT).

For example:
XS ttbar semileptonic = XS * 2( (BR(el+ve)+BR(µ+vµ)+BR(tau+vtau)) * BR(hadronic))

##### Measuring X-Sections of Samples on McM / DAS client #####
https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToGenXSecAnalyzer#Automated_scripts_to_compute_the

Find the MINIAOD dataset name on DAS client if you are using the command that requires the dataset.

Find your sample on McM. Find the MINIAOD prepID of the chain if you are using the command that requires the prepID.
e.g.
https://cms-pdmv.cern.ch/mcm/requests?prepid=B2G-RunIIFall17MiniAOD-00045&page=0&shown=127

Two possible commands.

For dataset input:
$>./calculateXSectionAndFilterEfficiency.sh -f datasets.txt -c Moriond17 -d MINIAODSIM -n 1000000

For McM prepID input:
$>./calculateXSectionAndFilterEfficiency.sh -f datasets_mcm.txt  -m -n 1000000


Once the file has finished running you will get output shown beneath:

'''
25-Jan-2018 10:39:51 CET  Closed file root://eoscms.cern.ch//eos/cms/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/50000/A20F16FB-E8EC-E711-9135-FA163E7ABC84.root

------------------------------------
GenXsecAnalyzer:
------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Overall cross-section summary
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Process		xsec_before [pb]		passed	nposw	nnegw	tried	nposw	nnegw 	xsec_match [pb]			accepted [%]	 event_eff [%]
0		5.458e+01 +/- 5.467e-03		107372	107298	74	340198	339942	256	1.723e+01 +/- 4.359e-02		31.6 +/- 0.1	31.6 +/- 0.1
1		2.528e+02 +/- 2.532e-02		368328	368014	314	1570294	1569026	1268	5.929e+01 +/- 8.582e-02		23.5 +/- 0.0	23.5 +/- 0.0
2		4.187e+02 +/- 4.192e-02		350533	350275	258	2602604	2600354	2250	5.641e+01 +/- 8.894e-02		13.5 +/- 0.0	13.5 +/- 0.0
3		3.610e+02 +/- 3.616e-02		174286	174135	151	2242464	2240173	2291	2.806e+01 +/- 6.473e-02		7.8 +/- 0.0	7.8 +/- 0.0
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Total		1.087e+03 +/- 6.112e-02		1000519	999722	797	6755560	6749495	6065	1.610e+02 +/- 1.491e-01		14.8 +/- 0.0	14.8 +/- 0.0
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Before matching: total cross section = 1.087e+03 +- 6.112e-02 pb
After matching: total cross section = 1.610e+02 +- 1.491e-01 pb
Matching efficiency = 0.1 +/- 0.0   [TO BE USED IN MCM]
Filter efficiency (taking into account weights)= (998925) / (998925) = 1.000e+00 +- 0.000e+00
Filter efficiency (event-level)= (1.00052e+06) / (1.00052e+06) = 1.000e+00 +- 0.000e+00    [TO BE USED IN MCM]

After filter: final cross section = 1.610e+02 +- 1.491e-01 pb
After filter: final fraction of events with negative weights = 7.966e-04 +- 1.386e-09
After filter: final equivalent lumi for 1M events (1/fb) = 6.190e+00 +- 3.548e-02
'''


Before matching = cross-section the output before any matching or filter.
After matching = cross-section after matching but before filter.
Filter efficiency = efficiency of any filter.
After filter = final cross-section.


##### Filter Efficiencies ####

https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideSubgroupMC#Measure_filter_efficiencies_and

You can use the script from the twiki

OR

Just find your sample on McM e.g.
https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_test/B2G-RunIIFall17DRPremix-00045

The find the associated wmLHEGS in the chain (I think you basically need the LHE fragments):
https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_test/B2G-RunIIFall17wmLHEGS-00071

Click on the circle with the tick on it to get the test command. Copy and paste into a shell script and this should run.
In the output look for the string "GenXsecAnalyzer".
Below this you have a detailed logging of the needed quantities.

##################################################################


#######################  Creation of Ntuples  #####################

######################## # Installation ######################## 
See base installation here https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagAnalyzer

######################## # Producing the trees locally ########################
  $> cmsRun runBTagAnalyzer_cfg.py defaults=Moriond18 runOnData=False miniAOD=True useTTbarFilter=True maxEvents=-1

- Will run locally the analyzer for testing purposes.

- Note the new way to initialize all the relevant defaults using defaults=Moriond18

- useTTbarFilter option adds ttbar variables: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagAnalyzer#Variables_produced_with_use_ttba

- Need to comment out line 10 from RecoBTag/Configuration/python/RecoBTag_cff.py to remove DeepFlavour.

- Need to update '/CMSSW_X_Y_Z/src/RecoBTag/PerformanceMeasurements/python/TTbarSelectionProducer_cfi.py’ and related cc and h files are where you set triggers, base selections etc.

- May be worth ensuring the default sample in PerformanceMeasurments/python/defaults/Moriond18.py to a larger sample otherwise bugs may arise in crab jobs instead of local test e.g.:
            - 'inputFiles' : ['/store/mc/RunIIFall17MiniAOD/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/60000/002E7FEA-16E0-E711-922D-0242AC130002.root'],

########################
# Producing trees on the grid
########################

Ensure you have a working CMSSW environment, you voms proxy set up and crab setup:
$> cd CMSSW_XXXXX/src/
$> cmsenv
$> voms-proxy-init --voms cms
$> source /cvmfs/cms.cern.ch/crab3/crab.sh

- Change into the directory from which you will submit the jobs e.g. ‘CMSSW_8_0_26_patch1/src/RecoBTag/PerformanceMeasurements/test/ttbar/‘
- Running the following commands will submit the samples described in the data/samples_XXX.json to run on the grid:

$> python submitToGrid.py -j data/samples_Run2016_8023.json -c ${CMSSW_BASE}/ src/RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py -l /afs/ cern.ch/work/f/fromeo/CMSSW_8_0_23_BTagAna/src/RecoBTag/ PerformanceMeasurements/test/ttbar/data/ Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt -s
$> python submitToGrid.py -j data/syst_samples_Run2016_80X_25ns.json -c $ {CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/ runBTagAnalyzer_cfg.py -s ```

- Partial submission can be made using the option:  -o csv_list
- In the samples list json, you should consider to leave the names of the MC and data samples as they were originally, with the prefix MC13TeV_ and Data13TeV_, to assure consistency toward all the steps of the measurement.

- If using new samples the cross-sections may need changing/adding. The best place to find these is in the McM. Look in the wmLHEGS production files. Note, these may have been produced with a different CMSSW - don’t worry.

- To process more quickly the MC jobs can edit number of units per job.
- Need to update the good runs list in submitToGrid.py : parser.add_option('-l', '--lumi', dest='lumiMask', help='json with list of good lumis', default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/ Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt')
- Adjust —lfn as needed. 
- Githash variable creates the folder where you save the crab output. You can change it into a more friendly name.
- githash=commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1]  can be e.g.


########################
# Copy trees
########################

$> source EosCpDir.sh //For MC samples source
$>EosCpDirMuEG.sh //For data samples 

UPDATE: For this you can use python script '/src/RecoBTag/PerformanceMeasurements/test/ttbar/scripts/xrdcp_script.py'
        This ensures the TTBar samples are split into training and testing at this stage so the pickle files for the
        luminosity scaling are calculated for the number of events in the split samples.

- This script substitute Pedros' checkProductionIntegrity.py script. It is easier (to me) and one can run as soon as ntuple production starts to end, to move from crab output directories to a more simple directory structure which can be which is easily parsed by the local analysis. Once production is finished you can delete original crab directories in EOS.
- Make sure that the number of files in the new simplified path is equal to the one produced with crab. Then you have the option to delete the crab folder note that, e.g. for data, the same datasetname may have from series: MuonEG may have datasample B,C.
- In this case create dedicated script to do the copy (see EosCpDirMuEG.sh). It is always a good practice to check that the files have been properly copied (especially for data).
- Use the script CheckNt.sh. If some files are missing, or new files are arriving from crab, you can use the script MissingFiles_EosCpDir.sh.
- Ensure that the size of the directories are the same for old and copied directory.
########################
# Pileup estimation
########################
 $> python runPileupEstimation.py --json data/Run2016_MuEGBC_processedLumis.json
or
$> python runPileupEstimation.py --json data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RUNB.txt


- Pileup weights correct for the differences between the modeled pileup scenario in simulation and what is measured in data.

- Number of primary vertices is the observable used to determine the weights.

- To create the pileup distribution in data runPileupEstimation.py uses the pileupCalc.py tool.

- The pileup distribution in MC is hard coded. Script uses mixing package to directly access MC pileup distributions but need to update line 6.
  MC: SimGeneral.MixingModule.mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi

- To calculate the final weight for an event:
      - Normalise both distributions.
      - Divide data by MC.
      - Eventually one determines the y-value of the resulting distribution @ the x-value corresponding to the true number of interactions of the MC event.
      - In ntuples:
                    nPUTrue c.f. PileupSummaryInfo object ipu->getTrueNumInteractions() function.
                    nPU c.f. PileupSummaryInfo object ipu->getPU_NumInteractions() function.
      - The number of pileup in any given event must be an integer.
      - nPUTrue is not an integer and is drawn from the full PU distribution input that is a poisson mean of the distribubtion one gets nPU from.
      - nPU is the actual number of interactions in the event.

- Script to produce a ROOT file under data with the pileup distributions and the weights for a conservative +/-10% variation of the central minBias xsec value assumed.

- /CMSSW_8_0_26_patch1/src/RecoBTag/PerformanceMeasurements/test/ttbar/data/Run2016_lumiSummary.json is an example of a json file with processed runs (--json option).

- Need to use crab to produce the json file with the lumi sections analysed:
    $> crab report dataset_task_name

- Need to update line 6.
  Data: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt

- If script crashes while using pileUpCalc.py because 'Run X not found in Lumi/Pileup input file.  Check your files!', just ensure that for the the --puJson file is the correct one for the era of data you are running.

- Needed to create seperated Cert<XXX>.json files from Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunBCDEF.txt file. Used das query to find the run numbers for the individual runs e.g.
  $>das_client.py --query="run dataset=/MuonEG/Run2017C-17Nov2017-v1/MINIAOD"
  Then just removed all other runs that werent in the list returned by the das query.

######################### Creation of Rtuples ########################
###########################
# Local Rtuples creation
###########################

IMPORTANT: Split your ttbar sample to ensure you don't overtrain your classifier.
           - For this you can use python script '/src/RecoBTag/PerformanceMeasurements/test/ttbar/scripts/xrdcp_script.py'


>>WARNING!!!! When you split your ttbar samples up for training/testing BDT purposes, one must delete and remake the pickle files.<<

$>python runTTbarAnalysis.py -i /store/group/phys_btag/Commissioning/TTbar/KinMethod_StructuredDir_08Jan2017/ -o output_ntuples/ -j data/samples_Run2017_ttbar_SL.json -n 50

-  It runs the local analysis to produce the Rtples for the efficiency measurement. MC will be weighted by cross section.
- Option -n indicates how many threads should be used.  
- The first time the script will produce a pickle file with the weights to be used according to the number of files found, xsec specified in the json file.
- It is advised to run a samples on a single thread the first time you run (suggest running just one event).

- runTTbarAnalysis.py creates a TTbarEventAnalysis object.
- Object defined in TTbarEventAnalysis.cc/.h
- Runs a few functions to set values for e.g triggers/mva/etc.
- Then runs sequence of functions: prepareOutput, processFile, finalizeOutput.
- Output stored in directory named after '-o' command line option.
- Input directory after '-i' option.
- Runs on samples in .json file.


- StoreTools.py:
  produceNormalizationCache() function loops over a list of samples and produces a cache file to normalize MC.  
  getEOSlslist() function takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended.
  It prints out s.th. like "Produced normalization cache (data/.xsecweights.pck)" the first time you run and the pickle file takes a while.
  It seems that when you prepare the pickle file for the first time it is important to start with one MC sample only at the beginning and do not use -n (but it will finish with seg fault, though should not affect the results) or use -n 2
  After the first MC sample is processed, you can use even higher -n (e.g. -n 8), but it seems you have to run the sample 1by1.

Things to check:
- Make sure that the number of rtplas equals the num of ntplas or use the script MissingRtplas.sh to adjust. Sometimes this is needed, because the processing of ntplas->rtplas is blocked Note

- In case you update the trees, xsec or lumi you have to remove by hand the pickle file. 
  rm nohup.out
   rm data/.xsecweights.pck 

- If output directory needs clearing or need to create a new one:
  eos rm -r /store/group/phys_btag/Commissioning/TTbar/XXXXX/<output_file_dir>
  xrd eoscms mkdir /store/group/phys_btag/Commissioning/TTbar/XXXXX/<output_file_dir>

- Modifications to storeTools.py :
  #eos_cmd = '/afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select'
  eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'  

- Modifications to runTTbarAnalysis.py
  Comment the section remove command from 'merge' section.
  Eos is mounted on lxplus therefore can be treated like a regular file system i.e. can use commands like hadd, rm mkdir ..... just so long as you have the correct path name. In this case '/eos/cms/' needed to be added in front of '/store/.....'.
  Added '/eos/cms/' to outFile name as path was not correct for writing to eos.

- Modifications to  TTbarEventAnalysis.h:
  Ensure latest and greatest object corrections e.g. /src/RecoBTag/PerformanceMeasurements/test/ttbar/data/Summer15_25nsV5_DATA_Uncertainty_AK4PFchs.txt
  Are you using the correct pileup reweighting file: look for 'jecUncUrl' variable?
  Ensure triggers are set up correctly.

- Modifications to TTbarEventAnalysis.cc: Trigger and Letpon SF
  TString outfile = "root://eoscms//eos/cms/"+outFile; outF_=TFile::Open(outfile,"RECREATE"); //outF_=TFile::Open(outFile,"RECREATE");
  Add modifications to selection:
	if(mll<90) continue;  (inline with analysis note)
	if(lp4[0].Pt()<25 || lp4[1].Pt()<25) continue; //According to trigger plateau
  Check functions are up-to-date:
    getTriggerEfficiency()
    getLeptonSelectionEfficiencyScaleFactor() - Lepton scale factors are applied on the event weight. The scale factors vary wrt lepton pt and eta.
    getJetResolutionScales() https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution

  Trigger Efficiencies:
  MuonEG triggers:
  Single Muons triggers:
  EGamma triggers:

  Lepton Selection Efficiencies:
  - This can be a bit of a ball ache.
  - For preliminary results you may have to use an the latest version of the scale factors even if they are a little old.
  - Ivan stated there should be preliminary values for 2017 data.

  Muons:
    - For muons they provide a .json: https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults.
    - Considering there is only one lepton scale factor here I believe this must be the product of isolation and ID scale factors.
    -

  Electrons:
    - Have a look here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors.
    - Considering there is only one lepton scale factor here I believe this must be the product of reconstruction and ID scale factors.
    - No .json provided by EGamma group so either ask another ttbar analyst if they have the SF's or can get from .root files.


########################
# Control plots
########################
 
$>python plotter.py -i output_ntuples_outTheBoxCode_MVA/ -j data/samples_Run2017_prelim_4plotter.json

 - Makes control plots and stores all in a ROOT file. Different options may be passed to filter plots, and show differently the plots.

- When merging rootplates, be careful because if different sample names are called similarly (eg tW and atW)  you can risk doing double-merging.

- When scaling to a given luminosity (--lumi) ensure that you enter the luminosity in inverse picobarns i.e. 1000 pb-1 not 1 fb-1.


=================================
      2017 Data Luminosity
=================================
$>export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH

$>brilcalc lumi -i data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RUNB.txt -u /pb
Luminosity: 684.765754486

$>brilcalc lumi -i data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RUNC.txt -u /pb
Luminosity: 9755.779650397

$>brilcalc lumi -i data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunD.txt -u /pb
Luminosity: 4319.797670448

$>brilcalc lumi -i data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RUNE.txt -u /pb
Luminosity: 9423.592284089

$>brilcalc lumi -i data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunF.txt -u /pb
Luminosity: 13568.206

$>brilcalc lumi -i data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunBCDEF.txt -u /pb
Luminosity: 37752.141522137

########################
# Train KinDiscriminator
########################

$> sh KIN_runClassifier.sh 

- After running the local analysis use the kin tree stored in the ttbar sample to train a kinematics discriminator for b-jets in ttbar events. The script compiles and runs KIN_trainClassifier.C which should be modified in case different trainings are required. 
- Can script to create parent directory (shouldnt be needed, added code to create dir if doesn't exist):
  $> source CreateKinDir.sh
- Rearrange new folders created by KIN_runClassifier.sh.
- Don't forget to split your ttbar sample to ensure you don't overtrain your classifier.

!!!!WARNING!!!!
- The number of events used in the training has been optimised to get the best separation between signal and background while maintaining a flat background.
- Changing this number can drastically change the shape of the kinematic discriminant.

########################
### Re-run Ntuples (includes new mva training in Rootplas)
########################
 $>  python runTTbarAnalysis.py -i /store/group/phys_btag/Commissioning/TTbar/ TTbarNt8023/ -o /store/group/phys_btag/Commissioning/TTbar/TTbarRt8023MVA/ -j data/samples_Run2016_8023.json --tmvaWgts data/KIN/ -n 4 

- Re-run the analysis to store the KIN discriminator value per jet and merge appropriately.
- Make sure that the number of rtplas equals the num of ntplas.
- Sometimes this is needed, because the processing of ntplas->rtplas is blocked


########################
# Fit
########################

>>> Prepare and Perform Fit <<<

$> python Templated_btagEffFitter.py --taggers data/taggers_Run2017.json --inDir output_ntuples_RunBCDEFpuw_addDeepCSV_2018-01-17_MVA/ --outDir fit_dir --recycleTemplates -n 500

Also has option to run channel:
--channels <one of the below identifiers>
- 121=ee,-143=emu,-169=mumu

IMPORTANT: Note some selection that was not present in TTBarEventAnalysis is present in Templated_btagEffFitter.py:
           e.g. if njets<2 or njets>4 : continue

- Note: Make sure that all the rootplas are correctly merged: check their size, to confirm that, because of similar names, some files contain duplicated rootplas.
- check that the names of the rootfile are properly saved, with the prefix MC13TeV_ and Data13TeV_.
- Make sure that the tagger json, e.g. taggers_Run2017.json, contains the right implementation of names and WP cuts.
- The fit is done in two parts:
  1.) Templated_btagEffFitter.py where the templates are prepared and the fit is required.
      Possible changes to be made:
	       Pt bins - 'jetpt': [(30,50),(50,70),(70,100),(100,140),(140,200),(200,300),(300,600)],
         Number of jets - if njets<2 or njets>4 : continue 
         flavour groups - flavourGroups=[ ['b'], ['c','other'] ]
  2.) TTbarSFbFitTools.cc where the Roofit implementation of fit.
      Possible changes to be made:
	       - Comment out and add the following:
           //RooRealVar *sfVar = new RooRealVar(Form("sf%d",i),title,1,0,2);
           //sfVars.add(*sfVar); 
           //if(i!=idxOfInterest) sfVar->setConstant(true); 	//Account for MisTag rate 
           RooRealVar *sfVar; 
           if(i==0) sfVar = new RooRealVar(Form("sf%d",i),title,1,0,2); 
           if(i==1) sfVar = new RooRealVar(Form("sf%d",i),title,mistagrate,0,2);
           sfVars.add(*sfVar); 
           if(i!=idxOfInterest) sfVar->setConstant(true);
        - label->DrawLatex(0.45,0.94,i==0 ? "CSVv2 M" : "CSVv2 M");
          Already edited to dynamically draw correct WP and jet pt range.
          Will need to add something here for deepcsv to ensure we get the plots.

    3.) TTbarSFbFitTools.h:
          -Float_t lumi=37.75214

- In particular the value of mistag rate has to be defined before, according to tagger, WP, and jetpt case. See dedicated lines.
- This means you also have to modify TTbarSFbFitTools.cc


>>> Post Fit <<<
- Performed by ' PostFitter.sh '.
- Creates a summary of plots of efficiency measurements and sf values printed in a tex file Note.
- It may be the script stops because of some latex reason. If so keep pressing return until the script finishes.
- Among the other results, a csv.root files is produced at some points.
- The histos contained in it are normalised to 1/pb.
- In particular b, c, other_pass0_slice0 contains the inclusive distributions, before any cuts on b tag discriminator WPs.
- So from b, c, other_pass0_slice0->Integral() you can get the relative contribution of b, c, other cases
- And from b, c, other_pass0_slice0->DrawNormalized() you can compare the templates for the b, c, other cases.
- In general we expect the b and c templates to be comparable.


>>> Prepare Results for Distribution <<<
- Performed by 'prepare_csv.sh '.
- Prepares the .csv file to share the results with other people Note
- You really need to read carefully all the instructions at the beginning of 'prepare_csv.sh'.
- There are a lot of small things to do to ensure the script performs.


>>> Systematic Errors <<<
- Use script:
$> root -l SommaQuadratura.cc
- SommaQuadratura.cc is macro that helps to
  a. evaluate the relative error of the measurement
  b. understanding how the relative error of the measurement evolves when adding new sys
  c. recalculate the values plus/minus the total errors for each pT bin, when new sys are added
  d. calculate the values of central measurement plus/minus the error due to a new systeamatic.
- This is useful when you have to use the systeamtic error of a discriminator (e.g. calculated for the csv) to another discriminator (e.g. cmva). Otherwise these values are taken directly from the measured values of the systematic.

- Also useful:
$>  ./prepare_OtherSys.sh 
- Takes as input the values of c. and d. of the previous script The values of c. have to substitute the current value in the .csv file The values of d. have to be added to it

Notes:
- I usually proceed like this 
1. Add sfvals and sferrs in SommaQuadratura.cc, as indicated in SommaQuadratura.cc and check the relative error, just to have a feeling that things are consistent with what we expect. This part is done where there is the line //Relative error 
2. Add the errors for the different new systematics in the dedicated array. Check how the Relative Error evolves when adding a systematics. See the line //Relative error evolution 
3. Take the new val in SommaQuadratura.cc (see line with //New vals) and add them to prepare_OtherSys.sh Add all the sys values taken from SommaQuadratura.cc, before running prepare_OtherSys.sh Note that you have to change the sign + and - for each variation up and down. I suggest you do toterr down and up, and then for all the other systematics you do all down first and all up after. This is how you naturally fill prepare_OtherSys.sh and I found it less confusing 4. Substitute the values of up and down in the default .csv file Add all the new down and all the new up sys in the default .csv file.


########################
#     Systematics
########################

>>> ttbar weight systematics <<<

Some useful reading:

[1]Evaluating Theoretical Uncertainties with EDM event weights and which weights are stored in MC:
https://twiki.cern.ch/twiki/bin/view/CMS/ScaleAndPDFUncertaintiesFromEventWeights

[2]Accessing the event weights from an EDM file with CMSSW:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW

[3]Top systematics:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Factorization_and_renormalizatio

[4]PSWeights:
https://github.com/cms-sw/cmssw/pull/21477

[5]Weight number convention
https://twiki.cern.ch/twiki/bin/view/CMS/TopModGen#Event_Generation




In this code, generator weights are stored in the output from the crab jobs (jobs run by the BTagAnalyzer).
The weights are stored in a ttree branch called 'ttbar_w' which is filled @ around line 1460:

CMSSW_X_Y_Z/src/RecoBTag/PerformanceMeasurements/plugins/BTagAnalyzer.cc

The way in which the weights are accessed and retrieved is explained in [2].








>>> PDF weights in my sample? <<<

To know which integer or string corresponds to which weight use the package here:

git clone https://github.com/kdlong/TheoreticalUncertainties

It needs to be used inside a CMSSW environment but otherwise runs pretty well out of the box:

python getWeightInfoFromEDMFile.py -f /store/mc/RunIIFall17MiniAOD/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/60000/002E7FEA-16E0-E711-922D-0242AC130002.root

The printout will take the form:

<weightgroup combine="envelope" name="scale_variation">
<weight id="1001"> muR=1 muF=1 </weight>
<weight id="1002"> muR=1 muF=2 </weight>
<weight id="1003"> muR=1 muF=0.5 </weight>
<weight id="1004"> muR=2 muF=1 </weight>

<weightgroup combine="hessian" name="PDF_variation">
<weight id="2001"> PDF set = 260001 </weight>
<weight id="2002"> PDF set = 260002 </weight>

So for example, the command above gave me:

<weightgroup combine="envelope" name="scale_variation">
<weight id="1001"> lhapdf=306000 renscfact=1d0 facscfact=1d0 </weight>
<weight id="1002"> lhapdf=306000 renscfact=1d0 facscfact=2d0 </weight>
<weight id="1003"> lhapdf=306000 renscfact=1d0 facscfact=0.5d0 </weight>

<weightgroup combine="hessian" name="PDF_variation1">
<weight id="2000"> lhapdf=306000 </weight>
<weight id="2001"> lhapdf=306001 </weight>
<weight id="2002"> lhapdf=306002 </weight>
<weight id="2003"> lhapdf=306003 </weight>
<weight id="2004"> lhapdf=306004 </weight>




>>> Systematic Samples <<<

Systematic samples:
https://cms-pdmv.cern.ch/mcm/requests?range=TOP-RunIIFall17wmLHEGS-00042,TOP-RunIIFall17wmLHEGS-00065&page=0&shown=127



#############################################################################################################################################################################################################
                                                                                FIN
#############################################################################################################################################################################################################


############################
#       brilCalc.py
############################

- Working with centrally installed version on lxplus.
- Calculated luminosities come from json files (good run list) that are output by crab jobs.

$>export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH

- Can then run python script using input .json file containing runs.

$>brilcalc lumi -i data/PileUp/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_RunD.txt


############################
#          GitHub
############################
Git basics:
https://github.com/GarageGames/Torque2D/wiki/Cloning-the-repo-and-working-with-Git

Red Alert Commands:
https://github.com/blog/2019-how-to-undo-almost-anything-with-git


>> Undo Local changes <<
Accidentally deleted your file and saved it/Disk mount has cleared all code from file (can happen!):
git checkout -- <bad filename>



############################
#         GENERAL
############################
USEFUL COMMAND for finding phrase in file in large directory e.g.:
$>grep nPUtrue */*
