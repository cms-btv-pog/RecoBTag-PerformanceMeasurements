import FWCore.ParameterSet.Config as cms

ttbarselectionfilter = cms.EDFilter('TTbarSelectionFilter',
                                    channel        = cms.InputTag("ttbarselectionproducer"),
                                    selectChannels = cms.vint32(-11*11,-13*13,-11*13),
                                    selectAll      = cms.bool(False)
                                    )
