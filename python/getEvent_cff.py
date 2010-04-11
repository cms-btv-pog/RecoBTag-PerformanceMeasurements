# Should monitor what is recommended in
# https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes
# For now used the old
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TRKPromptFeedBack#Event_and_track_selection_recipe
#

import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

###### BSC trigger bit

from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *

from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
bit40 = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('40 AND NOT (36 OR 37 OR 38 OR 39)'))

###### Physics declared trigger bit

from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
physDecl = hltHighLevelDev.clone(HLTPaths = ['HLT_PhysicsDeclared'], HLTPathsPrescales = [1])

###### One primary vertex

oneGoodVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

###### No scraping

noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

##########Counter for the number of processed events
eventCountProducer = cms.EDProducer("EventCountProducer")


getEventDATA = cms.Sequence(physDecl*bit40*eventCountProducer*noScraping*oneGoodVertexFilter)
getEventMC = cms.Sequence(eventCountProducer*noScraping*oneGoodVertexFilter)

