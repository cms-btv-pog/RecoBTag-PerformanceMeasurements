import FWCore.ParameterSet.Config as cms


OperatingPoints = cms.EDAnalyzer("OperatingPoints",

                    bTagCutList = cms.untracked.VPSet(
    cms.PSet(
	    collection = cms.untracked.InputTag('trackCountingHighEffBJetTags'),
	    alias = cms.untracked.string('TCHE'),
	    MinimumDiscriminator = cms.untracked.double(-1),
	    MaximumDiscriminator = cms.untracked.double(15),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.1),
		    name = cms.untracked.string('TCHEL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.01),
		    name = cms.untracked.string('TCHEM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.001),
		    name = cms.untracked.string('TCHET')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('trackCountingHighPurBJetTags'),
	    alias = cms.untracked.string('TCHP'),
	    MinimumDiscriminator = cms.untracked.double(-1),
	    MaximumDiscriminator = cms.untracked.double(15),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.1),
		    name = cms.untracked.string('TCHPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.01),
		    name = cms.untracked.string('TCHPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.001),
		    name = cms.untracked.string('TCHPT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('jetProbabilityBJetTags'),
	    alias = cms.untracked.string('JP'),
	    MinimumDiscriminator = cms.untracked.double(0),
	    MaximumDiscriminator = cms.untracked.double(1.5),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.1),
		    name = cms.untracked.string('JPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.01),
		    name = cms.untracked.string('JPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.001),
		    name = cms.untracked.string('JPT')
		    )
		) ),
        cms.PSet(
                collection = cms.untracked.InputTag('jetBProbabilityBJetTags'),
                            alias = cms.untracked.string('JBP'),
                            MinimumDiscriminator = cms.untracked.double(0),
                            MaximumDiscriminator = cms.untracked.double(1.5),
                            OperatingPoints = cms.untracked.VPSet(
    cms.PSet(
                        cut = cms.untracked.double(0.1),
                                            name = cms.untracked.string('JBPL')
                                            ),
                                    cms.PSet(
                        cut = cms.untracked.double(0.01),
                                            name = cms.untracked.string('JBPM')
                                            ),
                                    cms.PSet(
                        cut = cms.untracked.double(0.001),
                                            name = cms.untracked.string('JBPT')
                                            )
                                    ) ),
    
    cms.PSet(
                collection = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),
                alias = cms.untracked.string('SSV'),
                MinimumDiscriminator = cms.untracked.double(0),
                MaximumDiscriminator = cms.untracked.double(7),
                OperatingPoints = cms.untracked.VPSet(
    cms.PSet(
    cut = cms.untracked.double(0.1),
    name = cms.untracked.string('SSVL')
    ),
    cms.PSet(
    cut = cms.untracked.double(0.01),
    name = cms.untracked.string('SSVM')
    ),
    cms.PSet(
    cut = cms.untracked.double(0.001),
    name = cms.untracked.string('SSVT')
    )
    ) ),
    cms.PSet(
    collection = cms.untracked.InputTag('combinedSecondaryVertexBJetTags'),
                            alias = cms.untracked.string('CSV'),
                            MinimumDiscriminator = cms.untracked.double(0),
                            MaximumDiscriminator = cms.untracked.double(1),
                            OperatingPoints = cms.untracked.VPSet(
                    cms.PSet(
                        cut = cms.untracked.double(0.1),
                                            name = cms.untracked.string('CSVL')
                                            ),
                                    cms.PSet(
                        cut = cms.untracked.double(0.01),
                                            name = cms.untracked.string('CSVM')
                                            ),
                                    cms.PSet(
                        cut = cms.untracked.double(0.001),
                                            name = cms.untracked.string('CSVT')
                                            )
                                    ) ),
    
    
    ),
                    OperatingPointsbyMistagRate = cms.bool (True),
                    ApplyJetCorrections = cms.bool (True),
                    jetCorrectionsLabel = cms.string("L2L3JetCorrectorIC5Calo"),
                    jetcuts = cms.PSet(
	MaxEta = cms.double(2.5),
	MinPt = cms.double(30.0),
	MinNtracks = cms.int32(1),
	MinTrkPt = cms.double(1.),
	MinNjets = cms.int32(1)
	),
                    flavourSource = cms.InputTag('IC5byValAlgo'),
                    Jets = cms.string('iterativeCone5CaloJets'),
                    bTagTrackEventIPtagInfos = cms.string('impactParameterTagInfos'),
                    outputFile = cms.untracked.string('OperatingPoints.root'),
                    debug = cms.untracked.bool (False)

                    )


                    
