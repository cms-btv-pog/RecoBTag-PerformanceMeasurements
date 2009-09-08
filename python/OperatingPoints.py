import FWCore.ParameterSet.Config as cms


# To estimate Operating Points by udsg mistag rate
# below cut = mistag rate
EstimateByMistagRate = cms.untracked.PSet(

    OperatingPointsList = cms.untracked.VPSet(
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
	    ) )    
)
)

# Operating Points for 31X
# where cut is applied to default discriminator
OperatingPoints31X = cms.untracked.PSet(

    OperatingPointsList = cms.untracked.VPSet(
    cms.PSet(
	    collection = cms.untracked.InputTag('trackCountingHighEffBJetTags'),
	    alias = cms.untracked.string('TCHE'),
	    MinimumDiscriminator = cms.untracked.double(-1),
	    MaximumDiscriminator = cms.untracked.double(15),
	    OperatingPoints = cms.untracked.VPSet(
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
	    OperatingPoints = cms.untracked.VPSet(
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
	    OperatingPoints = cms.untracked.VPSet(
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
                            OperatingPoints = cms.untracked.VPSet(
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
                OperatingPoints = cms.untracked.VPSet(
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
	    OperatingPoints = cms.untracked.VPSet(
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
		) )

    )
)
