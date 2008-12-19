import FWCore.ParameterSet.Config as cms

caloJetCollectionClone = cms.EDProducer("CaloJetShallowCloneProducer",
    src = cms.InputTag("iterativeCone5CaloJets")
)

tagJet = cms.EDFilter("JetFlavourIdentifier",
    jets = cms.InputTag("caloJetCollectionClone"),
    debug = cms.bool(False),
    coneSizeToAssociate = cms.double(0.3),
    vetoFlavour = cms.vstring(),
    physicsDefinition = cms.bool(False)
)



