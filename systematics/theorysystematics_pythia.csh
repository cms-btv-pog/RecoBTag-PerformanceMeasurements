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
    set cfi=PythiaDefault.cfi
else
    set cfi=Pythia_${label}_${dir}.cfi
endif


### create cfi file:
echo "Creating file ${cfi}..."
if (-f "${cfi}") rm ${cfi}
cat > ${cfi} <<EOF
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
EOF

#more ${cfi}


echo " "
echo "Config file: ${cfi}"
#echo "Output file: ${output}"
#echo "Feed it into your SIM or FastSim job."
