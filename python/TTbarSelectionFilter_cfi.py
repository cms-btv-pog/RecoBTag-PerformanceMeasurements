import FWCore.ParameterSet.Config as cms

ttbarselectionfilter = cms.EDFilter('TTbarSelectionFilter',
                                    selectChannels = cms.vint32(-11*11,-13*13,-11*13), #dilepton channels
				    #selectChannels = cms.vint32(-11,11,-13,13), #semileptonic channels
                                    selectAll      = cms.bool(False)
                                    )
