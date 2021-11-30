import FWCore.ParameterSet.Config as cms

def customizePFPatLikeJets(process, runCalo=True, runPuppi=True, runPF=True):
    # set some default collection variables
    pfjets =                "hltAK4PFJets"                                  #original ak4PFJetsCHS
    puppijets =             "hltAK4PFPuppiJets"                                  #original ak4PFJetsCHS
    pfjetsCorrected =       "hltAK4PFJetsCorrected"                         #original ak4PFJetsCHS
    calojets =              "hltAK4CaloJets"                                #original ak4CaloJets
    calojetsCutted =        "hltSelectorCentralJets30L1FastJeta2p5"
    PFDeepCSVTags =         "hltDeepCombinedSecondaryVertexBPFPatJetTags"   #original pfDeepCSVJetTags
    PFPuppiDeepCSVTags =    "hltDeepCombinedSecondaryVertexBPFPuppiPatJetTags"   #original pfDeepCSVJetTags
    CaloDeepCSVTags =       "hltDeepCombinedSecondaryVertexCaloPatBJetTags"
    PFDeepFlavourTags =     "hltPFDeepFlavourPatJetTags"                       #original pfDeepFlavourJetTagsSlimmedDeepFlavour
    PFPuppiDeepFlavourTags ="hltPFPuppiDeepFlavourJetTags"                       #original pfDeepFlavourJetTagsSlimmedDeepFlavour
    rho =                   "hltFixedGridRhoFastjetAll"                     #original fixedGridRhoFastjetAll
    hltVertices =           "hltVerticesPFFilter"                           #original offlinePrimaryVertices
    siPixelClusters =       "hltSiPixelClusters"                            #original siPixelClusters
    ecalRecHit =            "hltEcalRecHit"                                 #original ecalRecHit
    hbhereco =              "hltHbhereco"                                   #original hbhereco
    hfreco =                "hltHfreco"                                     #original hfreco
    horeco =                "hltHoreco"                                     #original horeco
    rpcRecHits =            "hltRpcRecHits"                                 #original rpcRecHits
    # tracks =                "hltMergedTracks"                               #original generalTracks
    tracks =                "hltPFMuonMerging"                               #original generalTracks
    # tracks =                "hltPixelTracks"                               #original generalTracks
    payload =               "AK4PFHLT"                                      #original AK4PFchs
    payloadPuppi =          "AK4PFPuppiHLT"                                      #original AK4PFchs
    particleFlow =          "hltParticleFlow"                               #original particleFlow
    puppi =                 "hltPFPuppi"                                    #original puppi
    puppiNoLep =            "hltPFPuppiNoLep"                               #original puppiNoLep
    beamSpot =              "hltOnlineBeamSpot"                             #original offlineBeamSpot
    caloTower =             "hltTowerMakerForAll"

    # clone and modify the HLT BTV sequence/producers to remove the jet pt and eta selections from "jetsForBtag" and replace with pfjets
    process.hltDeepCombinedSecondaryVertexBPFPatJetTags = process.hltDeepCombinedSecondaryVertexBJetTagsPF.clone(
        src = cms.InputTag( "hltDeepCombinedSecondaryVertexBJetPatTagInfos" )
    )
    process.hltDeepCombinedSecondaryVertexBJetPatTagInfos = process.hltDeepCombinedSecondaryVertexBJetTagsInfos.clone(
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
        + process.hltDeepCombinedSecondaryVertexBJetPatTagInfos
        + process.hltDeepCombinedSecondaryVertexBPFPatJetTags
    )

    # do the same for PuppiJets

    process.hltDeepCombinedSecondaryVertexBPFPuppiPatJetTags = process.hltDeepCombinedSecondaryVertexBJetTagsPF.clone(
        src = cms.InputTag( "hltDeepCombinedSecondaryVertexBPuppiJetPatTagInfos" )
    )
    process.hltDeepCombinedSecondaryVertexBPuppiJetPatTagInfos = process.hltDeepCombinedSecondaryVertexBJetTagsInfos.clone(
        svTagInfos = cms.InputTag( "hltDeepSecondaryVertexPFPuppiPatTagInfos" )
    )
    process.hltDeepSecondaryVertexPFPuppiPatTagInfos = process.hltDeepSecondaryVertexTagInfosPF.clone(
        trackIPTagInfos = cms.InputTag( "hltDeepBLifetimePFPuppiPatTagInfos" ),
        weights = cms.InputTag(puppi)
    )
    process.hltDeepBLifetimePFPuppiPatTagInfos = process.hltDeepBLifetimeTagInfosPF.clone(
        jets = cms.InputTag( puppijets )
    )
    process.HLTBtagDeepCSVSequencePFPuppiPat = cms.Sequence(
        process.hltVerticesPF
        + process.hltVerticesPFSelector
        + process.hltVerticesPFFilter
        + process.hltDeepBLifetimePFPuppiPatTagInfos
        + process.hltDeepInclusiveVertexFinderPF
        + process.hltDeepInclusiveSecondaryVerticesPF
        + process.hltDeepTrackVertexArbitratorPF
        + process.hltDeepInclusiveMergedVerticesPF
        + process.hltDeepSecondaryVertexPFPuppiPatTagInfos
        + process.hltDeepCombinedSecondaryVertexBPuppiJetPatTagInfos
        + process.hltDeepCombinedSecondaryVertexBPFPuppiPatJetTags
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
            cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagInfos"),
            cms.InputTag("hltDeepSecondaryVertexPFPatTagInfos"),
            cms.InputTag("hltPFDeepFlavourPatTagInfos"),
        ),
        trackAssociationSource = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPF"),
    )
    process.hltPatJetsPuppi = patJets.clone(
        JetFlavourInfoSource = cms.InputTag("hltPatJetFlavourAssociationPuppi"),
        JetPartonMapSource = cms.InputTag("hltPatJetFlavourAssociationLegacyPuppi"),
        addJetID = cms.bool(False),
        addTagInfos = cms.bool(True),
        discriminatorSources = cms.VInputTag(
            cms.InputTag(PFPuppiDeepCSVTags,"probb"),cms.InputTag(PFPuppiDeepCSVTags,"probc"),cms.InputTag(PFPuppiDeepCSVTags,"probudsg"),
            # cms.InputTag(PFPuppiDeepCSVTags,"probbb"), # hltDeepCSV: probb = probb +probbb
            cms.InputTag(PFPuppiDeepFlavourTags,"probb"), cms.InputTag(PFPuppiDeepFlavourTags,"probc"), cms.InputTag(PFPuppiDeepFlavourTags,"probg"),
            cms.InputTag(PFPuppiDeepFlavourTags,"problepb"), cms.InputTag(PFPuppiDeepFlavourTags,"probbb"), cms.InputTag(PFPuppiDeepFlavourTags,"probuds"),
        ),
        embedGenPartonMatch = cms.bool(False),
        genJetMatch = cms.InputTag("hltPatJetGenJetMatchPuppi"),
        genPartonMatch = cms.InputTag("hltPatJetPartonMatchPuppi"),
        jetChargeSource = cms.InputTag("patJetPuppiCharge"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("hltPatJetCorrFactorsPuppi")),
        jetIDMap = cms.InputTag("hltAk4JetID"),
        jetSource = cms.InputTag(puppijets),
        tagInfoSources = cms.VInputTag(
            cms.InputTag("hltDeepBLifetimePFPuppiPatTagInfos"),
            cms.InputTag("hltDeepCombinedSecondaryVertexBPuppiJetPatTagInfos"),
            cms.InputTag("hltDeepSecondaryVertexPFPuppiPatTagInfos"),
            cms.InputTag("hltPFPuppiDeepFlavourTagInfos"),
        ),
        trackAssociationSource = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPFPuppi"),
    )
    process.hltPatJetsCalo = patJets.clone(
        JetFlavourInfoSource = cms.InputTag("hltPatJetFlavourAssociationCalo"),
        JetPartonMapSource = cms.InputTag("hltPatJetFlavourAssociationLegacyCalo"),
        addAssociatedTracks = cms.bool(True),
        addBTagInfo = cms.bool(True),
        addDiscriminators = cms.bool(True),
        addEfficiencies = cms.bool(False),
        embedCaloTowers = cms.bool(True),
        # addGenJetMatch = cms.bool(True),
        # addGenPartonMatch = cms.bool(True),
        addJetCharge = cms.bool(False),
        addJetCorrFactors = cms.bool(False),
        # addJetFlavourInfo = cms.bool(True),
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
        # jetChargeSource = cms.InputTag("hltPatJetCharge"),
        # jetCorrFactorsSource = cms.VInputTag(cms.InputTag("hltPatJetCorrFactors")),
        # jetIDMap = cms.InputTag("hltAk4JetID"),
        jetSource = cms.InputTag(calojetsCutted),
        tagInfoSources = cms.VInputTag(
            cms.InputTag("hltImpactParameterPatTagInfos"),
            cms.InputTag("hltDeepCombinedSecondaryVertexBJetCaloPatTagInfos"),
            cms.InputTag("hltInclusiveSecondaryVertexFinderPatTagInfos"),
            # cms.InputTag("hltImpactParameterTagInfos"),
            # cms.InputTag("hltInclusiveSecondaryVertexFinderTagInfos"),
            # cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsInfoCalo"),
            # cms.InputTag("hltDeepSecondaryVertexPFPuppiPatTagInfos"),
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
    process.hltPatJetFlavourAssociationPuppi = patJetFlavourAssociation.clone(
        bHadrons = cms.InputTag("hltPatJetPartons","bHadrons"),
        cHadrons = cms.InputTag("hltPatJetPartons","cHadrons"),
        jets = cms.InputTag(puppijets),
        leptons = cms.InputTag("hltPatJetPartons","leptons"),
        partons = cms.InputTag("hltPatJetPartons","physicsPartons"),
        weights = cms.InputTag(puppi)
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
    process.hltPatJetFlavourAssociationLegacyCalo = patJetFlavourAssociationLegacy.clone(
        srcByReference = cms.InputTag("hltPatJetPartonAssociationLegacyCalo")
    )
    process.hltPatJetFlavourAssociationLegacyPuppi = patJetFlavourAssociationLegacy.clone(
        srcByReference = cms.InputTag("hltPatJetPartonAssociationLegacyPuppi")
    )

    process.hltPatJetPartonAssociationLegacy = patJetPartonAssociationLegacy.clone(
        jets = cms.InputTag(pfjets),
        partons = cms.InputTag("hltPatJetPartonsLegacy")
    )
    process.hltPatJetPartonAssociationLegacyPuppi = patJetPartonAssociationLegacy.clone(
        jets = cms.InputTag(puppijets),
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
    process.hltPatJetGenJetMatchPuppi = patJetGenJetMatch.clone(
        matched = cms.InputTag("hltSlimmedGenJets"),
        src = cms.InputTag(puppijets)
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
    process.hltPatJetPartonMatchPuppi = patJetPartonMatch.clone(
        matched = cms.InputTag("hltPrunedGenParticles"),
        src = cms.InputTag(puppijets)
    )
    process.hltPatJetPartonMatchCalo = patJetPartonMatch.clone(
        matched = cms.InputTag("hltPrunedGenParticles"),
        src = cms.InputTag(calojetsCutted)
    )

    from PhysicsTools.PatAlgos.recoLayer0.jetTracksCharge_cff import patJetCharge
    process.hltPatJetCharge = patJetCharge.clone(
        src = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPF"),
    )
    process.patJetPuppiCharge = patJetCharge.clone(
        src = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPFPuppi"),
    )

    from RecoJets.JetAssociationProducers.ak4JTA_cff import ak4JetTracksAssociatorAtVertexPF
    process.hltAk4JetTracksAssociatorAtVertexPF = ak4JetTracksAssociatorAtVertexPF.clone(
        jets = cms.InputTag(pfjets),
        pvSrc = cms.InputTag(hltVertices),
        tracks = cms.InputTag(tracks),
    )
    process.hltAk4JetTracksAssociatorAtVertexPFPuppi = ak4JetTracksAssociatorAtVertexPF.clone(
        jets = cms.InputTag(puppijets),
        pvSrc = cms.InputTag(hltVertices),
        tracks = cms.InputTag(tracks),
    )
    process.hltAk4JetTracksAssociatorAtVertexCalo = ak4JetTracksAssociatorAtVertexPF.clone(
        jets = cms.InputTag(calojetsCutted),
        # pvSrc = cms.InputTag(hltVertices),
        pvSrc = cms.InputTag("hltVerticesL3"),
        # tracks = cms.InputTag(tracks),
        tracks = cms.InputTag("hltMergedTracksForBTag"),
    )

    from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi  import patJetCorrFactors
    process.hltPatJetCorrFactors = patJetCorrFactors.clone(
        payload = cms.string(payload),
        primaryVertices = cms.InputTag(hltVertices),
        rho = cms.InputTag(rho),
        src = cms.InputTag(pfjets),
    )
    process.hltPatJetCorrFactorsPuppi = patJetCorrFactors.clone(
        payload = cms.string(payloadPuppi),
        primaryVertices = cms.InputTag(hltVertices),
        rho = cms.InputTag(rho),
        src = cms.InputTag(puppijets),
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
    # from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags
    # process.hltPFDeepFlavourJetTags = pfDeepFlavourJetTags.clone(
    #     src = cms.InputTag("hltPFDeepFlavourTagInfos")
    # )
    # process.hltPFPuppiDeepFlavourJetTags = pfDeepFlavourJetTags.clone(
    #     src = cms.InputTag("hltPFPuppiDeepFlavourTagInfos")
    # )
    # from RecoBTag.FeatureTools.pfDeepFlavourTagInfos_cfi import pfDeepFlavourTagInfos
    # process.hltPFDeepFlavourTagInfos = pfDeepFlavourTagInfos.clone(
    #     candidates = cms.InputTag(particleFlow),
    #     jets = cms.InputTag(pfjets),
    #     fallback_puppi_weight = cms.bool(True),
    #     puppi_value_map = cms.InputTag(""),
    #     secondary_vertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
    #     shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagInfos"),
    #     vertex_associator = cms.InputTag("hltPrimaryVertexAssociation","original"),
    #     vertices = cms.InputTag(hltVertices)
    # )
    # process.hltPFPuppiDeepFlavourTagInfos = pfDeepFlavourTagInfos.clone(
    #     candidates = cms.InputTag(particleFlow),
    #     jets = cms.InputTag(puppijets),
    #     puppi_value_map = cms.InputTag(puppi),
    #     secondary_vertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
    #     shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBPuppiJetPatTagInfos"),
    #     vertex_associator = cms.InputTag("hltPrimaryVertexAssociationPuppi","original"),
    #     vertices = cms.InputTag(hltVertices)
    # )

    process.hltPrimaryVertexAssociationPat = process.hltPrimaryVertexAssociation.clone(
        jets = cms.InputTag(pfjets),
    )

    process.hltPFDeepFlavourPatTagInfos = process.hltPFDeepFlavourTagInfos.clone(
        jets = cms.InputTag(pfjets),
        shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagInfos"),
        vertex_associator = cms.InputTag("hltPrimaryVertexAssociationPat","original"),
    )

    process.hltPFDeepFlavourPatJetTags = process.hltPFDeepFlavourJetTags.clone(
        src = cms.InputTag("hltPFDeepFlavourPatTagInfos")
    )

    from RecoBTag.SecondaryVertex.candidateCombinedSecondaryVertexV2Computer_cfi import candidateCombinedSecondaryVertexV2Computer
    process.candidateCombinedSecondaryVertexV2Computer = candidateCombinedSecondaryVertexV2Computer.clone()

    # from PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi import primaryVertexAssociation
    # process.hltPrimaryVertexAssociation = primaryVertexAssociation.clone(
    #     jets = cms.InputTag(pfjets),
    #     particles = cms.InputTag(particleFlow),
    #     vertices = cms.InputTag(hltVertices),
    # )
    # process.hltPrimaryVertexAssociationPuppi = primaryVertexAssociation.clone(
    #     jets = cms.InputTag(puppijets),
    #     particles = cms.InputTag(particleFlow),
    #     vertices = cms.InputTag(hltVertices),
    # )
    process.HLTBtagDeepJetSequencePFPat = cms.Sequence(
        process.hltVerticesPF
        + process.hltVerticesPFSelector
        + process.hltVerticesPFFilter
        + process.hltDeepBLifetimePFPatTagInfos
        + process.hltDeepInclusiveVertexFinderPF
        + process.hltDeepInclusiveSecondaryVerticesPF
        + process.hltDeepTrackVertexArbitratorPF
        + process.hltDeepInclusiveMergedVerticesPF
        + process.hltDeepSecondaryVertexPFPatTagInfos
        + process.hltDeepCombinedSecondaryVertexBJetPatTagInfos
        + process.hltPrimaryVertexAssociationPat
        + process.hltPFDeepFlavourPatTagInfos
        + process.hltPFDeepFlavourPatJetTags
    )





    # create the final path
    if runPF:
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

            # *process.hltPrimaryVertexAssociationPat
            # *process.hltPFDeepFlavourPatTagInfos
            # *process.hltPFDeepFlavourPatJetTags
            * process.HLTBtagDeepJetSequencePFPat

            *process.hltPatJets
        )

    if runPuppi:
        process.MC_PuppiJetsMatchingPath = cms.Path(
            process.HLTAK4PFPuppiJetsSequence
            *process.HLTBtagDeepCSVSequencePFPuppiPat
            *process.hltPrunedGenParticlesWithStatusOne
            *process.hltPrunedGenParticles
            *process.hltPackedGenParticles
            *process.hltPatJetPartonMatchPuppi
            *process.hltSlimmedGenJets
            *process.hltAk4JetID
            *process.hltPatJetGenJetMatchPuppi
            *process.hltPatJetPartonsLegacy
            *process.hltPatJetPartonAssociationLegacyPuppi
            *process.hltPatJetFlavourAssociationLegacyPuppi
            *process.hltPatJetPartons
            *process.hltPatJetFlavourAssociationPuppi
            *process.hltAk4JetTracksAssociatorAtVertexPFPuppi
            *process.patJetPuppiCharge
            *process.hltPatJetCorrFactorsPuppi

            *process.hltPrimaryVertexAssociationPuppi
            *process.hltPFPuppiDeepFlavourTagInfos
            *process.hltPFPuppiDeepFlavourJetTags

            *process.hltPatJetsPuppi
        )

    if runCalo:
        process.MC_CaloJetsMatchingPath = cms.Path(
            process.HLTAK4CaloJetsCorrectionSequence
            *process.HLTBtagDeepCSVSequenceCaloPat
            *process.hltPrunedGenParticlesWithStatusOne
            *process.hltPrunedGenParticles
            *process.hltPackedGenParticles
            *process.hltPatJetPartonMatchCalo
            *process.hltPatJetGenJetMatchCalo
            *process.hltPatJetPartonsLegacy
            *process.hltPatJetPartons
            *process.hltSlimmedGenJets
            *process.hltPatJetPartonAssociationLegacyCalo
            *process.hltPatJetFlavourAssociationLegacyCalo
            *process.hltPatJetFlavourAssociationCalo
            *process.hltAk4JetTracksAssociatorAtVertexCalo

            *process.hltPatJetsCalo
        )

    if process.schedule_():
        if runPF: process.schedule.extend([process.MC_JetsMatchingPath])
        if runCalo: process.schedule.extend([process.MC_CaloJetsMatchingPath])
        if runPuppi: process.schedule.extend([process.MC_PuppiJetsMatchingPath])

    return process
