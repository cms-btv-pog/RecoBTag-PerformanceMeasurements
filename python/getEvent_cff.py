# Should monitor what is recommended in
# https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes
# 
#http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py?revision=1.15&view=markup
#

import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

###### BSC trigger bit

from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed

bit40_data = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True),
                                   L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')
                                   )

bit40_MC = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True),
                                 L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')
                                 )

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


getEventDATA = cms.Sequence(bit40_data*eventCountProducer*physDecl*noScraping*oneGoodVertexFilter)
getEventMC = cms.Sequence(bit40_MC*eventCountProducer*noScraping*oneGoodVertexFilter)

