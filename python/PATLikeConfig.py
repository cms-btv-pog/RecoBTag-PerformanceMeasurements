import FWCore.ParameterSet.Config as cms

def customizePFPatLikeJets(process, type = "AK4PFCHS"):
    # set some default collection variables
    pfjets =                "hltAK4PFJets"                                  #original ak4PFJetsCHS
    pfjetsCorrected =       "hltAK4PFJetsCorrected"                         #original ak4PFJetsCHS
    calojets =              "hltAK4CaloJets"                                #original ak4CaloJets
    calojetsCutted =        "hltSelectorCentralJets30L1FastJeta2p5"
    PFDeepCSVTags =         "hltDeepCombinedSecondaryVertexBPFPatJetTags"   #original pfDeepCSVJetTags
    CaloDeepCSVTags =       "hltDeepCombinedSecondaryVertexCaloPatBJetTags"
    PFDeepFlavourTags =     "hltPFDeepFlavourJetTags"                       #original pfDeepFlavourJetTagsSlimmedDeepFlavour
    rho =                   "hltFixedGridRhoFastjetAll"                     #original fixedGridRhoFastjetAll
    hltVertices =           "hltVerticesPFFilter"                           #original offlinePrimaryVertices
    siPixelClusters =       "hltSiPixelClusters"                            #original siPixelClusters
    ecalRecHit =            "hltEcalRecHit"                                 #original ecalRecHit
    hbhereco =              "hltHbhereco"                                   #original hbhereco
    hfreco =                "hltHfreco"                                     #original hfreco
    horeco =                "hltHoreco"                                     #original horeco
    rpcRecHits =            "hltRpcRecHits"                                 #original rpcRecHits
    tracks =                "hltMergedTracks"                               #original generalTracks
    payload =               "AK4PFHLT"                                      #original AK4PFchs
    particleFlow =          "hltParticleFlow"                               #original particleFlow
    puppi =                 "hltPFPuppi"                                    #original puppi
    puppiNoLep =            "hltPFPuppiNoLep"                               #original puppiNoLep
    beamSpot =              "hltOnlineBeamSpot"                             #original offlineBeamSpot

    # clone and modify the HLT BTV sequence/producers to remove the jet pt and eta selections from "jetsForBtag" and replace with pfjets
    process.hltDeepCombinedSecondaryVertexBPFPatJetTags = process.hltDeepCombinedSecondaryVertexBJetTagsPF.clone(
        src = cms.InputTag( "hltDeepCombinedSecondaryVertexBJetPatTagsInfos" )
    )
    process.hltDeepCombinedSecondaryVertexBJetPatTagsInfos = process.hltDeepCombinedSecondaryVertexBJetTagsInfos.clone(
        svTagInfos = cms.InputTag( "hltDeepSecondaryVertexPFPatTagInfos" )
    )
    process.hltDeepSecondaryVertexPFPatTagInfos = process.hltDeepSecondaryVertexTagInfosPF.clone(
        trackIPTagInfos = cms.InputTag( "hltDeepBLifetimePFPatTagInfos" )
    )
    process.hltDeepBLifetimePFPatTagInfos = process.hltDeepBLifetimeTagInfosPF.clone(
        jets = cms.InputTag( pfjets )
    )
    process.HLTBtagDeepCSVSequencePFPat = cms.Sequence(
        process.hltVerticesPF
        + process.hltVerticesPFSelector
        + process.hltVerticesPFFilter
        + process.hltDeepBLifetimePFPatTagInfos
        + process.hltDeepInclusiveVertexFinderPF
        + process.hltDeepInclusiveSecondaryVerticesPF
        + process.hltDeepTrackVertexArbitratorPF
        + process.hltDeepInclusiveMergedVerticesPF
        + process.hltDeepSecondaryVertexPFPatTagInfos
        + process.hltDeepCombinedSecondaryVertexBJetPatTagsInfos
        + process.hltDeepCombinedSecondaryVertexBPFPatJetTags
    )

    # do the same for caloJets
    process.hltDeepCombinedSecondaryVertexCaloPatBJetTags = process.hltDeepCombinedSecondaryVertexBJetTagsCalo.clone(
        src = cms.InputTag("hltDeepCombinedSecondaryVertexBJetCaloPatTagInfos"),
    )

    process.hltDeepCombinedSecondaryVertexBJetCaloPatTagInfos = process.hltDeepCombinedSecondaryVertexBJetTagsInfosCalo.clone(
        svTagInfos = cms.InputTag("hltInclusiveSecondaryVertexFinderPatTagInfos")
    )
    process.hltInclusiveSecondaryVertexFinderPatTagInfos = process.hltInclusiveSecondaryVertexFinderTagInfos.clone(
        trackIPTagInfos = cms.InputTag("hltImpactParameterPatTagInfos"),
    )

    process.hltImpactParameterPatTagInfos = process.hltImpactParameterTagInfos.clone(
        jetTracks = cms.InputTag("hltFastPixelBLifetimeL3AssociatorPat"),
    )

    process.hltSelectorCentralJets30L1FastJeta2p5 = process.hltSelectorCentralJets30L1FastJeta.clone(
        etaMax = cms.double(2.5),
        etaMin = cms.double(-2.5),
        src = cms.InputTag("hltSelectorJets20L1FastJet")
    )

    process.hltSelectorJets20L1FastJet = process.hltSelectorJets30L1FastJet.clone(
        etMin = cms.double(20.0),
    )

    process.hltFastPixelBLifetimeL3AssociatorPat = process.hltFastPixelBLifetimeL3Associator.clone(
        jets = cms.InputTag(calojetsCutted),
    )

    process.HLTBtagDeepCSVSequenceCaloPat = cms.Sequence(
        process.hltSelectorJets20L1FastJet
        +process.hltSelectorCentralJets30L1FastJeta2p5
        +process.HLTTrackReconstructionForBTag
        +process.hltVerticesL3
        +process.hltFastPixelBLifetimeL3AssociatorPat
        +process.hltImpactParameterPatTagInfos
        +process.hltInclusiveVertexFinder
        +process.hltInclusiveSecondaryVertices
        +process.hltTrackVertexArbitrator
        +process.hltInclusiveMergedVertices
        +process.hltInclusiveSecondaryVertexFinderPatTagInfos
        +process.hltDeepCombinedSecondaryVertexBJetCaloPatTagInfos
        +process.hltDeepCombinedSecondaryVertexCaloPatBJetTags
    )


    # create patJets  for ak4pfchs and all necessary missing inputs
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
    process.hltPatJets = patJets.clone(
        JetFlavourInfoSource = cms.InputTag("hltPatJetFlavourAssociation"),
        JetPartonMapSource = cms.InputTag("hltPatJetFlavourAssociationLegacy"),
        addJetID = cms.bool(False),
        addTagInfos = cms.bool(True),
        discriminatorSources = cms.VInputTag(
            cms.InputTag(PFDeepCSVTags,"probb"),cms.InputTag(PFDeepCSVTags,"probc"),cms.InputTag(PFDeepCSVTags,"probudsg"),
            # cms.InputTag(PFDeepCSVTags,"probbb"), # hltDeepCSV: probb = probb +probbb
            cms.InputTag(PFDeepFlavourTags,"probb"), cms.InputTag(PFDeepFlavourTags,"probc"), cms.InputTag(PFDeepFlavourTags,"probg"),
            cms.InputTag(PFDeepFlavourTags,"problepb"), cms.InputTag(PFDeepFlavourTags,"probbb"), cms.InputTag(PFDeepFlavourTags,"probuds"),
        ),
        embedGenPartonMatch = cms.bool(False),
        genJetMatch = cms.InputTag("hltPatJetGenJetMatch"),
        genPartonMatch = cms.InputTag("hltPatJetPartonMatch"),
        jetChargeSource = cms.InputTag("hltPatJetCharge"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("hltPatJetCorrFactors")),
        jetIDMap = cms.InputTag("hltAk4JetID"),
        jetSource = cms.InputTag(pfjets),
        tagInfoSources = cms.VInputTag(
            cms.InputTag("hltDeepBLifetimePFPatTagInfos"),
            cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagsInfos"),
            cms.InputTag("hltDeepSecondaryVertexPFPatTagInfos"),
            cms.InputTag("hltPFDeepFlavourTagInfos"),
        ),
        trackAssociationSource = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPF"),
    )
    process.hltPatJetsCalo = patJets.clone(
        JetFlavourInfoSource = cms.InputTag("hltPatJetFlavourAssociationCalo"),
        JetPartonMapSource = cms.InputTag("hltPatJetFlavourAssociationLegacy"),
            addAssociatedTracks = cms.bool(True),
            addBTagInfo = cms.bool(True),
            addDiscriminators = cms.bool(True),
            addEfficiencies = cms.bool(False),
            addGenJetMatch = cms.bool(True),
            addGenPartonMatch = cms.bool(True),
            addJetCharge = cms.bool(False),
            addJetCorrFactors = cms.bool(False),
            addJetFlavourInfo = cms.bool(True),
            addPartonJetMatch = cms.bool(False),
        addJetID = cms.bool(False),
        addTagInfos = cms.bool(True),
        discriminatorSources = cms.VInputTag(
            cms.InputTag(CaloDeepCSVTags,"probb"),cms.InputTag(CaloDeepCSVTags,"probc"),cms.InputTag(CaloDeepCSVTags,"probudsg"),
            # # cms.InputTag(PFDeepCSVTags,"probbb"), # hltDeepCSV: probb = probb +probbb
        ),
        embedGenPartonMatch = cms.bool(False),
        genJetMatch = cms.InputTag("hltPatJetGenJetMatchCalo"),
        genPartonMatch = cms.InputTag("hltPatJetPartonMatchCalo"),
        jetChargeSource = cms.InputTag("hltPatJetCharge"),
        # jetCorrFactorsSource = cms.VInputTag(cms.InputTag("hltPatJetCorrFactors")),
        # jetIDMap = cms.InputTag("hltAk4JetID"),
        jetSource = cms.InputTag(calojetsCutted),
        tagInfoSources = cms.VInputTag(
            cms.InputTag("hltDeepCombinedSecondaryVertexBJetCaloPatTagInfos"),
            cms.InputTag("hltImpactParameterPatTagInfos"),
            cms.InputTag("hltImpactParameterPatTagInfos"),
        ),
        trackAssociationSource = cms.InputTag("hltAk4JetTracksAssociatorAtVertexCalo"),
    )

    # for patJets
    from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import patJetFlavourAssociation,patJetPartons,patJetFlavourAssociationLegacy,patJetPartonAssociationLegacy,patJetPartonsLegacy
    process.hltPatJetFlavourAssociation = patJetFlavourAssociation.clone(
        bHadrons = cms.InputTag("hltPatJetPartons","bHadrons"),
        cHadrons = cms.InputTag("hltPatJetPartons","cHadrons"),
        jets = cms.InputTag(pfjets),
        leptons = cms.InputTag("hltPatJetPartons","leptons"),
        partons = cms.InputTag("hltPatJetPartons","physicsPartons"),
    )
    process.hltPatJetFlavourAssociationCalo = patJetFlavourAssociation.clone(
        bHadrons = cms.InputTag("hltPatJetPartons","bHadrons"),
        cHadrons = cms.InputTag("hltPatJetPartons","cHadrons"),
        jets = cms.InputTag(calojetsCutted),
        leptons = cms.InputTag("hltPatJetPartons","leptons"),
        partons = cms.InputTag("hltPatJetPartons","physicsPartons"),
    )
    process.hltPatJetPartons = patJetPartons.clone()

    process.hltPatJetFlavourAssociationLegacy = patJetFlavourAssociationLegacy.clone(
        srcByReference = cms.InputTag("hltPatJetPartonAssociationLegacy")
    )

    process.hltPatJetPartonAssociationLegacy = patJetPartonAssociationLegacy.clone(
        jets = cms.InputTag(pfjets),
        partons = cms.InputTag("hltPatJetPartonsLegacy")
    )
    process.hltPatJetPartonAssociationLegacyCalo = patJetPartonAssociationLegacy.clone(
        jets = cms.InputTag(calojetsCutted),
        partons = cms.InputTag("hltPatJetPartonsLegacy")
    )

    process.hltPatJetPartonsLegacy = patJetPartonsLegacy.clone(
        src = cms.InputTag("genParticles"),
    )

    from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetGenJetMatch
    process.hltPatJetGenJetMatch = patJetGenJetMatch.clone(
        matched = cms.InputTag("hltSlimmedGenJets"),
        src = cms.InputTag(pfjets)
    )
    process.hltPatJetGenJetMatchCalo = patJetGenJetMatch.clone(
        matched = cms.InputTag("hltSlimmedGenJets"),
        src = cms.InputTag(calojetsCutted)
    )

    from PhysicsTools.PatAlgos.slimming.slimmedGenJets_cfi import slimmedGenJets
    process.hltSlimmedGenJets = slimmedGenJets.clone(
        packedGenParticles = cms.InputTag("hltPackedGenParticles"),
        src = cms.InputTag("ak4GenJetsNoNu")
    )

    from PhysicsTools.PatAlgos.slimming.packedGenParticles_cfi import packedGenParticles
    process.hltPackedGenParticles = packedGenParticles.clone(
        inputCollection = cms.InputTag("hltPrunedGenParticlesWithStatusOne"),
        inputOriginal = cms.InputTag("genParticles"),
        map = cms.InputTag("hltPrunedGenParticles"),
    )

    from PhysicsTools.PatAlgos.slimming.genParticles_cff import prunedGenParticlesWithStatusOne
    from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import prunedGenParticles
    process.hltPrunedGenParticlesWithStatusOne = prunedGenParticlesWithStatusOne.clone(
        src = cms.InputTag("genParticles")
    )

    process.hltPrunedGenParticles = prunedGenParticles.clone(
        src = cms.InputTag("hltPrunedGenParticlesWithStatusOne")
    )

    from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetPartonMatch
    process.hltPatJetPartonMatch = patJetPartonMatch.clone(
        matched = cms.InputTag("hltPrunedGenParticles"),
        src = cms.InputTag(pfjets)
    )
    process.hltPatJetPartonMatchCalo = patJetPartonMatch.clone(
        matched = cms.InputTag("hltPrunedGenParticles"),
        src = cms.InputTag(calojetsCutted)
    )

    from PhysicsTools.PatAlgos.recoLayer0.jetTracksCharge_cff import patJetCharge
    process.hltPatJetCharge = patJetCharge.clone(
        src = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPF"),
    )

    from RecoJets.JetAssociationProducers.ak4JTA_cff import ak4JetTracksAssociatorAtVertexPF
    process.hltAk4JetTracksAssociatorAtVertexPF = ak4JetTracksAssociatorAtVertexPF.clone(
        jets = cms.InputTag(pfjets),
        pvSrc = cms.InputTag(hltVertices),
        tracks = cms.InputTag(tracks),
    )
    process.hltAk4JetTracksAssociatorAtVertexCalo = ak4JetTracksAssociatorAtVertexPF.clone(
        jets = cms.InputTag(calojetsCutted),
        pvSrc = cms.InputTag(hltVertices),
        tracks = cms.InputTag(tracks),
    )

    from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi  import patJetCorrFactors
    process.hltPatJetCorrFactors = patJetCorrFactors.clone(
        payload = cms.string(payload),
        primaryVertices = cms.InputTag(hltVertices),
        rho = cms.InputTag(rho),
        src = cms.InputTag(pfjets),
    )

    from RecoJets.JetProducers.ak4JetID_cfi import ak4JetID
    process.hltAk4JetID = ak4JetID.clone(
        ebRecHitsColl = cms.InputTag(ecalRecHit,"EcalRecHitsEB"),
        eeRecHitsColl = cms.InputTag(ecalRecHit,"EcalRecHitsEE"),
        hbheRecHitsColl = cms.InputTag(hbhereco),
        hfRecHitsColl = cms.InputTag(hfreco),
        hoRecHitsColl = cms.InputTag(horeco),
        rpcRecHits = cms.InputTag(rpcRecHits),
        src = cms.InputTag(calojets),
    )



    #### TAGGERS
    # run DeepFlavour for HLT
    from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags
    process.hltPFDeepFlavourJetTags = pfDeepFlavourJetTags.clone(
        src = cms.InputTag("hltPFDeepFlavourTagInfos")
    )
    from RecoBTag.FeatureTools.pfDeepFlavourTagInfos_cfi import pfDeepFlavourTagInfos
    process.hltPFDeepFlavourTagInfos = pfDeepFlavourTagInfos.clone(
        candidates = cms.InputTag(particleFlow),
        jets = cms.InputTag(pfjets),
        puppi_value_map = cms.InputTag(puppi),
        secondary_vertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
        shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagsInfos"),
        vertex_associator = cms.InputTag("hltPrimaryVertexAssociation","original"),
        vertices = cms.InputTag(hltVertices)
    )

    from RecoBTag.SecondaryVertex.candidateCombinedSecondaryVertexV2Computer_cfi import candidateCombinedSecondaryVertexV2Computer
    process.candidateCombinedSecondaryVertexV2Computer = candidateCombinedSecondaryVertexV2Computer.clone()

    from PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi import primaryVertexAssociation
    process.hltPrimaryVertexAssociation = primaryVertexAssociation.clone(
        jets = cms.InputTag(pfjets),
        particles = cms.InputTag(particleFlow),
        vertices = cms.InputTag(hltVertices),
    )

    from RecoParticleFlow.PFProducer.chargedHadronPFTrackIsolation_cfi import chargedHadronPFTrackIsolation
    process.hltChargedHadronPFTrackIsolation = chargedHadronPFTrackIsolation.clone(
        src = cms.InputTag(particleFlow)
    )

    # create the final path
    process.MC_JetsMatchingPath = cms.Path(
        process.HLTAK4PFJetsSequence
        *process.HLTBtagDeepCSVSequencePFPat
        *process.hltPrunedGenParticlesWithStatusOne
        *process.hltPrunedGenParticles
        *process.hltPackedGenParticles
        *process.hltPatJetPartonMatch
        *process.hltSlimmedGenJets
        *process.hltAk4JetID
        *process.hltPatJetGenJetMatch
        *process.hltPatJetPartonsLegacy
        *process.hltPatJetPartonAssociationLegacy
        *process.hltPatJetFlavourAssociationLegacy
        *process.hltPatJetPartons
        *process.hltPatJetFlavourAssociation
        *process.hltAk4JetTracksAssociatorAtVertexPF
        *process.hltPatJetCharge
        *process.hltPatJetCorrFactors

        *process.hltPrimaryVertexAssociation
        *process.hltChargedHadronPFTrackIsolation
        *process.hltPFDeepFlavourTagInfos
        *process.hltPFDeepFlavourJetTags

        *process.hltPatJets
        )

    process.MC_CaloJetsMatchingPath = cms.Path(
        process.HLTAK4CaloJetsCorrectionSequence
        *process.HLTBtagDeepCSVSequenceCaloPat
        *process.hltPrunedGenParticlesWithStatusOne
        *process.hltPrunedGenParticles
        *process.hltPackedGenParticles
        *process.hltPatJetPartonMatchCalo
        *process.hltPatJetGenJetMatchCalo
        *process.hltPatJetPartonsLegacy
        *process.hltSlimmedGenJets
        *process.hltPatJetPartonAssociationLegacyCalo
        *process.hltPatJetFlavourAssociationLegacy
        *process.hltPatJetFlavourAssociationCalo
        *process.hltAk4JetTracksAssociatorAtVertexCalo

        *process.hltPatJetsCalo
    )

    if process.schedule_():
        process.schedule.extend([process.MC_JetsMatchingPath])
        process.schedule.extend([process.MC_CaloJetsMatchingPath])

    return process
