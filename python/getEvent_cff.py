# Should monitor what is recommended in
# https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes
# 
#http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py?revision=1.15&view=markup
#

import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

### Remove events with anomalous HCAL noise ###
#process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *

###### One primary vertex

oneGoodVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

oneGoodVertexFilterMC = cms.EDFilter("VertexSelector",
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


getEventDATA = cms.Sequence(
    eventCountProducer*
    HBHENoiseFilter*
    noScraping*oneGoodVertexFilter)

getEventMC = cms.Sequence(
    eventCountProducer*noScraping*oneGoodVertexFilterMC)

