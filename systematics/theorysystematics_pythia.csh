#! /bin/csh

echo "This script will generate .cfi file containing the dedicated PYTHIA settings"

### default settings:

set parp61 = 0.25
set parp72 = 0.25
set parj81 = 0.25
set mstp3 = 1 #?

set parp67 = 1 #?
set parp71 = 4 #?

set parj42 = 0.52
set parj21 = 0.40

set parj54 = -0.031
set parj55 = -0.0041

set parp82 = 2.9

	
### menu:
echo "-------- Available systematics: ---------"
echo "0: default settings"
echo "1: lambda_qcd"
echo "2: Q2max"
echo "3: light quark fragmentation"
echo "4: heavy quark fragmentation"
echo "5: underlying event"
echo "---------------- Choose: ----------------"

set syst=$<

switch ( ${syst} )
    case 0:
	echo "you chose the default settings"
	goto filename
    case 1:
	echo "you chose lambdaqcd"
	set label = lambdaqcd
	breaksw
    case 2:
	echo "you chose Q2max"
	set label = q2max
	breaksw
    case 3:
	echo "you chose light quark fragmentation"
	set label = fragl
	breaksw
    case 4:
	echo "you chose heavy quark fragmentation"
	set label = fragh
	breaksw
    case 5:
	echo "you chose underlying event"
	set label = ue
	breaksw
    default:
	echo "Sorry, unrecognized option."
	exit
endsw


updown:
echo "---------------- Choose up/down: ----------------"
set dir=$<
switch ( ${dir} )
    case up:
	breaksw
    case down:
	breaksw
    default:
	echo "Sorry, unrecognized option. You must choose up or down."
	goto updown
endsw

### apply new settings:
switch ( ${syst} )
    case 1:
	if (${dir} == up) then
	    set parp61 = 0.35
	    set parp72 = 0.35
	    set parj81 = 0.35
	else
	    set parp61 = 0.15
	    set parp72 = 0.15
	    set parj81 = 0.15
	endif 
	breaksw
    case 2:
	if (${dir} == up) then
	    set parp67 = 4
	    set parp71 = 16
	else
	    set parp67 = 0.25
	    set parp71 = 1
	endif 
	breaksw
    case 3:
	if (${dir} == up) then
	    set parj42 = 0.56
	    set parj21 = 0.43
	else
	    set parj42 = 0.48
	    set parj21 = 0.37
	endif 
	breaksw
    case 4:
	if (${dir} == up) then
	    set parj54 = -0.020
	    set parj55 = -0.0037
	else
	    set parj54 = -0.042
	    set parj55 = -0.0045
	endif 
	breaksw
    case 5:
	if (${dir} == up) then
	    set parp82 = 3.4
	else
	    set parp82 = 2.4
	endif 
	breaksw
endsw

### name of the cfi which has to replace "Configuration/Generator/data/PythiaUESettings.cfi"
filename:
if (${syst} == 0) then
    set cfi=ttbar_fromScratch_fastsim_PythiaDefault.cfg
else
    set cfi=ttbar_fromScratch_fastsim_Pythia_${label}_${dir}.cfg
endif


### create cfi file:
echo "Creating file ${cfi}..."
if (-f "${cfi}") rm ${cfi}
cat > ${cfi} <<EOF
process PROD  = 
{

source = PythiaSource {
    untracked int32 pythiaPylistVerbosity = 0
    untracked bool pythiaHepMCVerbosity = false
    untracked int32 maxEventsToPrint = 0
    
    PSet PythiaParameters = {
	
	# This is a vector of ParameterSet names to be read, in this order
	vstring parameterSets = {
	    "pythiaUESettings", 
	    "processParameters"
	}
	
	vstring pythiaUESettings = {      
      'PARP(61)=${parp61}',
      'PARP(72)=${parp72}',
      'PARJ(81)=${parj81}',
      'MSTP(3)=${mstp3}',
      'PARP(67)=${parp67}',#
      'PARP(71)=${parp71}',
      'PARJ(42)=${parj42}',
      'PARJ(21)=${parj21}',
      'PARJ(54)=${parj54}',
      'PARJ(55)=${parj55}',
      'PARP(82)=${parp82}',#
	    'MSTJ(11)=3     ! Choice of the fragmentation function',
	    'MSTJ(22)=2     ! Decay those unstable particles',
	    'PARJ(71)=10 .  ! for which ctau  10 mm',
	    'MSTP(2)=1      ! which order running alphaS',
	    'MSTP(33)=0     ! no K factors in hard cross sections',
	    'MSTP(51)=7     ! structure function chosen',
	    'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default',
	    'MSTP(82)=4     ! Defines the multi-parton model',
	    'MSTU(21)=1     ! Check on possible errors during program execution',
	    'PARP(82)=1.9409   ! pt cutoff for multiparton interactions',
	    'PARP(89)=1960. ! sqrts for which PARP82 is set',
	    'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter',
	    'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter',
	    'PARP(90)=0.16  ! Multiple interactions: rescaling power',
	    'PARP(67)=2.5    ! amount of initial-state radiation',
	    'PARP(85)=1.0  ! gluon prod. mechanism in MI',
	    'PARP(86)=1.0  ! gluon prod. mechanism in MI',
	    'PARP(62)=1.25   ! ',
	    'PARP(64)=0.2    ! ',
	    'MSTP(91)=1     !',
	    'PARP(91)=2.1   ! kt distribution',
	    'PARP(93)=15.0  ! '
	}
	
	vstring processParameters = {
	    'MSEL=0                  ! User defined processes',
	    'MSUB(81) = 1            ! qqbar to QQbar',
	    'MSUB(82) = 1            ! gg to QQbar',
	    'MSTP(7) = 6             ! flavour = top',
	    'PMAS(6,1) = 175.        ! top quark mass'     
	}
    }
}


    # The number of events to be processed.
    untracked PSet maxEvents = {untracked int32 input = 830000}
    
    service =  RandomNumberGeneratorService {
	# This is to initialize the random engine of the source
	untracked uint32 sourceSeed = 123456789
	# This is to initialize the random engines of Famos
	PSet moduleSeeds =
	{
	    untracked uint32 VtxSmeared = 123456789
	    untracked uint32 famosPileUp = 918273
	    untracked uint32 famosSimHits = 13579
	    untracked uint32 siTrackerGaussianSmearingRecHits = 24680
	    untracked uint32 caloRecHits = 654321
	    untracked uint32 paramMuons = 54525
	}
    }
    
    // If you want to use the (CMS default) CLHEP random, 
    // set these ones to false
    replace famosPileUp.UseTRandomEngine = true
    replace famosSimHits.UseTRandomEngine = true
    replace siTrackerGaussianSmearingRecHits.UseTRandomEngine = true
    replace caloRecHits.UseTRandomEngine = true
    replace paramMuons.UseTRandomEngine = true

   
    # Famos sequences
    include "FastSimulation/Configuration/data/FamosSequences.cff"
    // If you want to turn on/off pile-up
    replace famosPileUp.PileUpSimulator.averageNumber = 0.0
    // Parametrized magnetic field
    replace VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = true
    // You may not want to simulate everything for your study
    replace famosSimHits.SimulateCalorimetry = true
    replace famosSimHits.SimulateTracking = true


    # Tracker MisAlignement 
    include "FastSimulation/Configuration/data/MisAlignment.cff" 
    replace misalignedTrackerGeometry.applyAlignment = true
    replace trackerAlignment.toGet = {
	{ string record = "TrackerAlignmentRcd" string tag = "Tracker10pbScenario150" },
	{ string record = "TrackerAlignmentErrorRcd" string tag = "Tracker10pbScenarioErrors150" }
    }

    # ECAL miscalibration. 
    # include "FastSimulation/Configuration/data/MisCalibration.cff"
	 
    ##replacement for ECAL miscalibration	
    include "CondCore/DBCommon/data/CondDBSetup.cfi"

    es_source ecalConditions = PoolDBESSource  
    { 
        using CondDBSetup
        VPSet toGet = {
	    { string record = "EcalIntercalibConstantsRcd"
	        string tag = "EcalIntercalibConstants_miscalibcsa07"}
        }
        string connect = "frontier://cms_conditions_data/CMS_COND_16X_ECAL"
        untracked bool siteLocalConfig = true
        string timetype = "runnumber"
    }

    replace caloRecHits.RecHitsFactory.doMiscalib=true
    replace ecalConditions.catalog="relationalcatalog_frontier://cms_conditions_data/CMS_COND_16X_FRONTIER"
	
    # HCAL miscalibration 
    # 1) RMS (1.0 means 10% RMS miscalibration, 0.5 means 5%, 2.0 means 20%)
    # Default is 0.0 (i.e., no miscalibration)
    replace hcalRecHits.Refactor = 1.0
    # 2) Bias (1.0 means no bias, 1.1 means 10% positive bias)
    # Default is 1.0 (i.e., no bias)
    replace hcalRecHits.Refactor_mean = 0.95	 	


    service = Timing { }
    
    path p1 = { 
	famosWithEverything
    }

    # To write out events (not need: FastSimulation _is_ fast!)
    include "FastSimulation/Configuration/data/EventContent.cff"
    module o1 = PoolOutputModule { 
	using AODSIMEventContent
	untracked string fileName = "ttbar_default.root" 
    }
    endpath outpath = { o1 }
    
    # Keep the logging output to a nice level #
    include "FWCore/MessageService/data/MessageLogger.cfi"
//    replace MessageLogger.destinations = {"detailedInfo.txt"}
    
}		

EOF

#more ${cfi}


echo " "
echo "Config file: ${cfi}"
#echo "Output file: ${output}"
#echo "Feed it into your SIM or FastSim job."
