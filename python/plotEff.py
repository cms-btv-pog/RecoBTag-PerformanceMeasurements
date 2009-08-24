import FWCore.ParameterSet.Config as cms


plotEff = cms.EDAnalyzer("plotEff",

                    bTagCutList = cms.untracked.VPSet(
    cms.PSet(
	    collection = cms.untracked.InputTag('trackCountingHighEffBJetTags'),
	    alias = cms.untracked.string('TCHE'),
	    MinimumDiscriminator = cms.untracked.double(-1),
	    MaximumDiscriminator = cms.untracked.double(15),
	    plotEff = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(1.87),
		    name = cms.untracked.string('TCHEL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(3.72),
		    name = cms.untracked.string('TCHEM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(14.2),
		    name = cms.untracked.string('TCHET')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('trackCountingHighPurBJetTags'),
	    alias = cms.untracked.string('TCHP'),
	    MinimumDiscriminator = cms.untracked.double(-1),
	    MaximumDiscriminator = cms.untracked.double(15),
	    plotEff = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(1.33),
		    name = cms.untracked.string('TCHPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(2.13),
		    name = cms.untracked.string('TCHPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(3.35),
		    name = cms.untracked.string('TCHPT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('jetProbabilityBJetTags'),
	    alias = cms.untracked.string('JP'),
	    MinimumDiscriminator = cms.untracked.double(0),
	    MaximumDiscriminator = cms.untracked.double(1.5),
	    plotEff = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.228),
		    name = cms.untracked.string('JPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.483),
		    name = cms.untracked.string('JPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.685),
		    name = cms.untracked.string('JPT')
		    )
		) ),
        cms.PSet(
                collection = cms.untracked.InputTag('jetBProbabilityBJetTags'),
                            alias = cms.untracked.string('JBP'),
                            MinimumDiscriminator = cms.untracked.double(0),
                            MaximumDiscriminator = cms.untracked.double(1.5),
                            plotEff = cms.untracked.VPSet(
    cms.PSet(
                        cut = cms.untracked.double(1.15),
                                            name = cms.untracked.string('JBPL')
                                            ),
                                    cms.PSet(
                        cut = cms.untracked.double(1.9),
                                            name = cms.untracked.string('JBPM')
                                            ),
                                    cms.PSet(
                        cut = cms.untracked.double(2),
                                            name = cms.untracked.string('JBPT')
                                            )
                                    ) ),
    
    cms.PSet(
                collection = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),
                alias = cms.untracked.string('SSV'),
                MinimumDiscriminator = cms.untracked.double(0),
                MaximumDiscriminator = cms.untracked.double(7),
                plotEff = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(1.22),
		    name = cms.untracked.string('SSVL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(1.8),
		    name = cms.untracked.string('SSVM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(3.53),
		    name = cms.untracked.string('SSVT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('combinedSecondaryVertexBJetTags'),
	    alias = cms.untracked.string('CSV'),
	    MinimumDiscriminator = cms.untracked.double(0),
	    MaximumDiscriminator = cms.untracked.double(1),
	    plotEff = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.422),
		    name = cms.untracked.string('CSVL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.815),
		    name = cms.untracked.string('CSVM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.944),
		    name = cms.untracked.string('CSVT')
		    )
		) ),
    
    
    ),
                    plotEffbyMistagRate = cms.bool (True),
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
                    outputFile = cms.untracked.string('plotEff.root'),
                    debug = cms.untracked.bool (False)

                    )


                    
