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
		cut = cms.untracked.double(0.02),
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
    cms.PSet(
	collection = cms.untracked.InputTag('softMuonBJetTags'),
	alias = cms.untracked.string('SMT'),
	MinimumDiscriminator = cms.untracked.double(0),
	MaximumDiscriminator = cms.untracked.double(1),
	OperatingPoints = cms.untracked.VPSet(
	    cms.PSet(
		cut = cms.untracked.double(0.012),
		name = cms.untracked.string('SMTL')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.01),
		name = cms.untracked.string('SMTM')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.003),
		name = cms.untracked.string('SMTT')
		)
	    ) ),
    cms.PSet(
	collection = cms.untracked.InputTag('softMuonByIP3dBJetTags'),
	alias = cms.untracked.string('SMTbyIP'),
	MinimumDiscriminator = cms.untracked.double(-10),
	MaximumDiscriminator = cms.untracked.double(30),
	OperatingPoints = cms.untracked.VPSet(
	    cms.PSet(
		cut = cms.untracked.double(0.012),
		name = cms.untracked.string('SMTbyIPL')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.01),
		name = cms.untracked.string('SMTbyIPM')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.003),
		name = cms.untracked.string('SMTbyIPT')
		)
	    ) ),
    cms.PSet(
	collection = cms.untracked.InputTag('softMuonByPtBJetTags'),
	alias = cms.untracked.string('SMTbyPt'),
	MinimumDiscriminator = cms.untracked.double(0),
	MaximumDiscriminator = cms.untracked.double(8),
	OperatingPoints = cms.untracked.VPSet(
	    cms.PSet(
		cut = cms.untracked.double(0.012),
		name = cms.untracked.string('SMTbyPtL')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.01),
		name = cms.untracked.string('SMTbyPtM')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.003),
		name = cms.untracked.string('SMTbyPtT')
		)
	    ) ),
    cms.PSet(
	collection = cms.untracked.InputTag('softElectronByIP3dBJetTags'),
	alias = cms.untracked.string('SETbyIP'),
	MinimumDiscriminator = cms.untracked.double(-10),
	MaximumDiscriminator = cms.untracked.double(30),
	OperatingPoints = cms.untracked.VPSet(
	    cms.PSet(
		cut = cms.untracked.double(0.012),
		name = cms.untracked.string('SETbyIPL')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.01),
		name = cms.untracked.string('SETbyIPM')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.003),
		name = cms.untracked.string('SETbyIPT')
		)
	    ) ),
    cms.PSet(
	collection = cms.untracked.InputTag('softElectronByPtBJetTags'),
	alias = cms.untracked.string('SETbyPt'),
	MinimumDiscriminator = cms.untracked.double(0),
	MaximumDiscriminator = cms.untracked.double(8),
	OperatingPoints = cms.untracked.VPSet(
	    cms.PSet(
		cut = cms.untracked.double(0.012),
		name = cms.untracked.string('SETbyPtL')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.01),
		name = cms.untracked.string('SETbyPtM')
		),
	    cms.PSet(
		cut = cms.untracked.double(0.003),
		name = cms.untracked.string('SETbyPtT')
		)
	    ) )
)
)

####################################################################################
# Operating Points for 31X
# where cut is applied to default discriminator
#
# operating points are estimated applying a loose taggability requirement,
# for 31X, we use the following samples:
# InclusiveMu5_Pt50 for the soft muon taggers
# QCD_Pt30 for the rest of taggers
#
OperatingPoints31X = cms.untracked.PSet(

    OperatingPointsList = cms.untracked.VPSet(
    cms.PSet(
	    collection = cms.untracked.InputTag('trackCountingHighEffBJetTags'),
	    alias = cms.untracked.string('TCHE'),
	    MinimumDiscriminator = cms.untracked.double(-1),
	    MaximumDiscriminator = cms.untracked.double(15),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(1.66),
		    name = cms.untracked.string('TCHEL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(2.96),
		    name = cms.untracked.string('TCHEM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(8.22),
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
		    cut = cms.untracked.double(1.13),
		    name = cms.untracked.string('TCHPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(1.85),
		    name = cms.untracked.string('TCHPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(3.08),
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
		    cut = cms.untracked.double(0.21),
		    name = cms.untracked.string('JPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.46),
		    name = cms.untracked.string('JPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.60),
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
		    cut = cms.untracked.double(0.96),
		    name = cms.untracked.string('JBPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(1.91),
		    name = cms.untracked.string('JBPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(2.06),
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
		    cut = cms.untracked.double(1.7),
		    name = cms.untracked.string('SSVM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(2.98),
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
		    cut = cms.untracked.double(0.37),
		    name = cms.untracked.string('CSVL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.74),
		    name = cms.untracked.string('CSVM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.92),
		    name = cms.untracked.string('CSVT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('softMuonBJetTags'),
	    alias = cms.untracked.string('SMT'),
	    MinimumDiscriminator = cms.untracked.double(0),
	    MaximumDiscriminator = cms.untracked.double(1),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.33),
		    name = cms.untracked.string('SMTL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.34),
		    name = cms.untracked.string('SMTM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.41),
		    name = cms.untracked.string('SMTT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('softMuonByIP3dBJetTags'),
	    alias = cms.untracked.string('SMTbyIP'),
	    MinimumDiscriminator = cms.untracked.double(-10),
	    MaximumDiscriminator = cms.untracked.double(30),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(8.22),
		    name = cms.untracked.string('SMTbyIPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(11.2),
		    name = cms.untracked.string('SMTbyIPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(32.2),
		    name = cms.untracked.string('SMTbyIPT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('softMuonByPtBJetTags'),
	    alias = cms.untracked.string('SMTbyPt'),
	    MinimumDiscriminator = cms.untracked.double(0),
	    MaximumDiscriminator = cms.untracked.double(8),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(1.96),
		    name = cms.untracked.string('SMTbyPtL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(2.09),
		    name = cms.untracked.string('SMTbyPtM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(3.42),
		    name = cms.untracked.string('SMTbyPtT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('softElectronByIP3dBJetTags'),
	    alias = cms.untracked.string('SETbyIP'),
	    MinimumDiscriminator = cms.untracked.double(-10),
	    MaximumDiscriminator = cms.untracked.double(30),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.74),
		    name = cms.untracked.string('SETbyIPL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(1.04),
		    name = cms.untracked.string('SETbyIPM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(3.41),
		    name = cms.untracked.string('SETbyIPT')
		    )
		) ),
    cms.PSet(
	    collection = cms.untracked.InputTag('softElectronByPtBJetTags'),
	    alias = cms.untracked.string('SETbyPt'),
	    MinimumDiscriminator = cms.untracked.double(0),
	    MaximumDiscriminator = cms.untracked.double(8),
	    OperatingPoints = cms.untracked.VPSet(
		cms.PSet(
		    cut = cms.untracked.double(0.66),
		    name = cms.untracked.string('SETbyPtL')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(0.73),
		    name = cms.untracked.string('SETbyPtM')
		    ),
		cms.PSet(
		    cut = cms.untracked.double(1.18),
		    name = cms.untracked.string('SETbyPtT')
		    )
		) )
    )

)
