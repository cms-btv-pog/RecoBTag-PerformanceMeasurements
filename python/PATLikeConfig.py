import FWCore.ParameterSet.Config as cms

def customizePFPatLikeJets(process):

    pfjets = "hltAK4PFJets" #original ak4PFJetsCHS
    pfjetsCorrected = "hltAK4PFJetsCorrected" #original ak4PFJetsCHS
    calojets = "hltAK4CaloJets" #original ak4CaloJets
    # PFDeepCSVTags = "hltDeepCombinedSecondaryVertexBJetTagsPF" # original: pfDeepCSVJetTags
    PFDeepCSVTags = "hltDeepCombinedSecondaryVertexBPFPatJetTags" # original: pfDeepCSVJetTags
    PFDeepFlavourTags = "hltPFDeepFlavourJetTags" # original: pfDeepFlavourJetTagsSlimmedDeepFlavour
    rho = "hltFixedGridRhoFastjetAll" #original fixedGridRhoFastjetAll
    hltVertices = "hltVerticesPFFilter" #original offlinePrimaryVertices
    # hltVerticesSlimmed = "hltVerticesPFFilter" #original offlineSlimmedPrimaryVertices
    hltVerticesSlimmed = "hltOfflineSlimmedPrimaryVertices" #original offlineSlimmedPrimaryVertices
    siPixelClusters = "hltSiPixelClusters" #original siPixelClusters
    ecalRecHit = "hltEcalRecHit" #original ecalRecHit
    hbhereco = "hltHbhereco" #original hbhereco
    hfreco = "hltHfreco" #original hfreco
    horeco = "hltHoreco" #original horeco
    rpcRecHits = "hltRpcRecHits" #original rpcRecHits
    tracks = "hltMergedTracks" #original generalTracks
    payload = "AK4PFHLT" #original AK4PFchs
    particleFlow = "hltParticleFlow" #original particleFlow
    puppi = "hltPFPuppi" #original puppi
    puppiNoLep = "hltPFPuppiNoLep" #original puppiNoLep
    beamSpot = "hltOnlineBeamSpot" #original offlineBeamSpot


    # clone and modify the HLT BTV sequence to remove the jet selection
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
        # jets = cms.InputTag( pfjetsCorrected )
        jets = cms.InputTag( pfjets )
    )
    process.HLTBtagDeepCSVSequencePFPat = cms.Sequence(
        process.hltVerticesPF
        + process.hltVerticesPFSelector
        + process.hltVerticesPFFilter
        # + process.hltPFJetForBtagSelector
        # + process.hltPFJetForBtag
        + process.hltDeepBLifetimePFPatTagInfos
        + process.hltDeepInclusiveVertexFinderPF
        + process.hltDeepInclusiveSecondaryVerticesPF
        + process.hltDeepTrackVertexArbitratorPF
        + process.hltDeepInclusiveMergedVerticesPF
        + process.hltDeepSecondaryVertexPFPatTagInfos
        + process.hltDeepCombinedSecondaryVertexBJetPatTagsInfos
        + process.hltDeepCombinedSecondaryVertexBPFPatJetTags
    )

    # process.hltSlimmedJets = cms.EDFilter("PATJetSelector",
    #     cut = cms.string('pt > 10'),
    #     cutLoose = cms.string(''),
    #     nLoose = cms.uint32(0),
    #     # src = cms.InputTag("hltUpdatedPatJetsTransientCorrectedSlimmedDeepFlavour")
    #     src = cms.InputTag("hltPatJets")
    # )
    # process.hltUpdatedPatJetsTransientCorrectedSlimmedDeepFlavour = cms.EDProducer("PATJetUpdater",
    #     # addBTagInfo = cms.bool(True),
    #     addBTagInfo = cms.bool(False),
    #     # addDiscriminators = cms.bool(True),
    #     addDiscriminators = cms.bool(False),
    #     addJetCorrFactors = cms.bool(True),
    #     addTagInfos = cms.bool(True),
    #     discriminatorSources = cms.VInputTag(
    #         # cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probb"), cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probc"), cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probg"),
    #         # cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","problepb"), cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probbb"), cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probuds"),
    #     ),
    #     jetCorrFactorsSource = cms.VInputTag(cms.InputTag("hltPatJetCorrFactorsTransientCorrectedSlimmedDeepFlavour")),
    #     # jetSource = cms.InputTag("hltUpdatedPatJetsSlimmedDeepFlavour"),
    #     jetSource = cms.InputTag("hltPatJets"),
    #     printWarning = cms.bool(True),
    #     tagInfoSources = cms.VInputTag("hltPixelClusterTagInfos"),
    # )

    # process.hltUpdatedPatJetsSlimmedDeepFlavour = cms.EDProducer("PATJetUpdater",
    #     # addBTagInfo = cms.bool(True),
    #     addBTagInfo = cms.bool(False),
    #     addDiscriminators = cms.bool(True),
    #     addJetCorrFactors = cms.bool(True),
    #     addTagInfos = cms.bool(False),
    #     # addTagInfos = cms.bool(True),
    #     jetCorrFactorsSource = cms.VInputTag(cms.InputTag("hltPatJetCorrFactorsSlimmedDeepFlavour")),
    #     # jetSource = cms.InputTag("hltSlimmedJetsNoDeepFlavour"),
    #     # jetSource = cms.InputTag("hltSelectedPatJets"),
    #     jetSource = cms.InputTag("hltPatJets"),
    #     printWarning = cms.bool(False),
    # )

    # process.hltPatJetCorrFactorsSlimmedDeepFlavour = cms.EDProducer("JetCorrFactorsProducer",
    #     emf = cms.bool(False),
    #     extraJPTOffset = cms.string('L1FastJet'),
    #     flavorType = cms.string('J'),
    #     levels = cms.vstring(),
    #     payload = cms.string(payload),
    #     primaryVertices = cms.InputTag(hltVerticesSlimmed),
    #     rho = cms.InputTag(rho),
    #     # src = cms.InputTag("hltSlimmedJetsNoDeepFlavour"),
    #     # src = cms.InputTag("hltSelectedPatJets"),
    #     src = cms.InputTag("hltPatJets"),
    #     useNPV = cms.bool(True),
    #     useRho = cms.bool(True)
    # )

    # process.hltSlimmedJetsNoDeepFlavour = cms.EDProducer("PATJetSlimmer",
    #     dropDaughters = cms.string('0'),
    #     dropJetVars = cms.string('1'),
    #     dropSpecific = cms.string('0'),
    #     dropTagInfos = cms.string('0'),
    #     dropTrackRefs = cms.string('1'),
    #     mixedDaughters = cms.bool(False),
    #     modifyJets = cms.bool(True),
    #     packedPFCandidates = cms.InputTag("hltPackedPFCandidates"),
    #     rekeyDaughters = cms.string('1'),
    #     # src = cms.InputTag("hltSelectedPatJets"),
    #     src = cms.InputTag("hltPatJets"),
    #     modifierConfig = cms.PSet(
    #         modifications = cms.VPSet()
    #     ),
    # )


    process.hltPackedPFCandidates = cms.EDProducer("PATPackedCandidateProducer",
        PuppiNoLepSrc = cms.InputTag(puppiNoLep),
        PuppiSrc = cms.InputTag(puppi),
        chargedHadronIsolation = cms.InputTag("hltChargedHadronPFTrackIsolation"),
        covariancePackingSchemas = cms.vint32(8, 264, 520, 776, 0),
        covarianceVersion = cms.int32(1),
        inputCollection = cms.InputTag(particleFlow),
        inputVertices = cms.InputTag(hltVerticesSlimmed),
        minPtForChargedHadronProperties = cms.double(3.0),
        minPtForTrackProperties = cms.double(0.95),
        originalTracks = cms.InputTag(tracks),
        originalVertices = cms.InputTag(hltVertices),
        pfCandidateTypesForHcalDepth = cms.vint32(),
        # secondaryVerticesForWhiteList = cms.VInputTag(cms.InputTag("inclusiveCandidateSecondaryVertices"), cms.InputTag("inclusiveCandidateSecondaryVerticesCvsL"), cms.InputTag("generalV0Candidates","Kshort"), cms.InputTag("generalV0Candidates","Lambda")),
        # secondaryVerticesForWhiteList = cms.VInputTag(cms.InputTag("inclusiveCandidateSecondaryVertices")),
        secondaryVerticesForWhiteList = cms.VInputTag(cms.InputTag("hltDeepInclusiveSecondaryVerticesPF")),
        storeHcalDepthEndcapOnly = cms.bool(False),
        storeTiming = cms.bool(False),
        timeMap = cms.InputTag(""),
        timeMapErr = cms.InputTag(""),
        vertexAssociator = cms.InputTag("hltPrimaryVertexAssociation","original")
    )

    # process.hltSelectedPatJets = cms.EDFilter("PATJetSelector",
    #     cut = cms.string('pt > 10'),
    #     cutLoose = cms.string(''),
    #     nLoose = cms.uint32(0),
    #     src = cms.InputTag("hltPatJets")
    # )

    process.hltPatJets = cms.EDProducer("PATJetProducer",
        JetFlavourInfoSource = cms.InputTag("hltPatJetFlavourAssociation"),
        JetPartonMapSource = cms.InputTag("hltPatJetFlavourAssociationLegacy"),
        addAssociatedTracks = cms.bool(True),
        addBTagInfo = cms.bool(True),
        addDiscriminators = cms.bool(True),
        addEfficiencies = cms.bool(False),
        addGenJetMatch = cms.bool(True),
        addGenPartonMatch = cms.bool(True),
        addJetCharge = cms.bool(True),
        addJetCorrFactors = cms.bool(True),
        addJetFlavourInfo = cms.bool(True),
        addJetID = cms.bool(False),
        # addPartonJetMatch = cms.bool(False),
        addPartonJetMatch = cms.bool(True),
        addResolutions = cms.bool(False),
        addTagInfos = cms.bool(True),
        discriminatorSources = cms.VInputTag(
            cms.InputTag(PFDeepCSVTags,"probb"),cms.InputTag(PFDeepCSVTags,"probc"),cms.InputTag(PFDeepCSVTags,"probudsg"),
            # cms.InputTag(PFDeepCSVTags,"probbb")
            cms.InputTag(PFDeepFlavourTags,"probb"), cms.InputTag(PFDeepFlavourTags,"probc"), cms.InputTag(PFDeepFlavourTags,"probg"),
            cms.InputTag(PFDeepFlavourTags,"problepb"), cms.InputTag(PFDeepFlavourTags,"probbb"), cms.InputTag(PFDeepFlavourTags,"probuds"),
        ),
        efficiencies = cms.PSet(),
        embedGenJetMatch = cms.bool(True),
        embedGenPartonMatch = cms.bool(False),
        embedPFCandidates = cms.bool(False),
        genJetMatch = cms.InputTag("hltPatJetGenJetMatch"),
        genPartonMatch = cms.InputTag("hltPatJetPartonMatch"),
        getJetMCFlavour = cms.bool(True),
        jetChargeSource = cms.InputTag("hltPatJetCharge"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("hltPatJetCorrFactors")),
        jetIDMap = cms.InputTag("hltAk4JetID"),
        jetSource = cms.InputTag(pfjets),
        partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
        resolutions = cms.PSet(),
        tagInfoSources = cms.VInputTag(
            cms.InputTag("hltPixelClusterTagInfos"),
            cms.InputTag("hltDeepBLifetimePFPatTagInfos"),
            cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagsInfos"),
            cms.InputTag("hltDeepSecondaryVertexPFPatTagInfos"),
            cms.InputTag("hltPFDeepFlavourTagInfos"),

        ),
        trackAssociationSource = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPF"),
        useLegacyJetMCFlavour = cms.bool(False),
        # userData = cms.PSet(
        #     userCands = cms.PSet(
        #         src = cms.VInputTag("")
        #     ),
        #     userClasses = cms.PSet(
        #         src = cms.VInputTag("")
        #     ),
        #     userFloats = cms.PSet(
        #         src = cms.VInputTag(
        #             cms.InputTag("pileupJetId","fullDiscriminant"), "QGTagger:qgLikelihood", "hfJetShowerShape:sigmaEtaEta", "hfJetShowerShape:sigmaPhiPhi", "caloJetMap:pt",
        #             "caloJetMap:emEnergyFraction"
        #         )
        #     ),
        #     userFunctionLabels = cms.vstring(),
        #     userFunctions = cms.vstring(),
        #     userInts = cms.PSet(
        #         src = cms.VInputTag(cms.InputTag("pileupJetId","fullId"), "hfJetShowerShape:centralEtaStripSize", "hfJetShowerShape:adjacentEtaStripsSize")
        #     )
        # )
    )

    # for patJets
    process.hltPatJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
        bHadrons = cms.InputTag("hltPatJetPartons","bHadrons"),
        cHadrons = cms.InputTag("hltPatJetPartons","cHadrons"),
        ghostRescaling = cms.double(1e-18),
        hadronFlavourHasPriority = cms.bool(False),
        jetAlgorithm = cms.string('AntiKt'),
        jets = cms.InputTag(pfjets),
        leptons = cms.InputTag("hltPatJetPartons","leptons"),
        partons = cms.InputTag("hltPatJetPartons","physicsPartons"),
        rParam = cms.double(0.4)
    )
    process.hltPatJetPartons = cms.EDProducer("HadronAndPartonSelector",
        fullChainPhysPartons = cms.bool(True),
        particles = cms.InputTag("genParticles"),
        partonMode = cms.string('Auto'),
        src = cms.InputTag("generator")
    )
    process.hltPatJetFlavourAssociationLegacy = cms.EDProducer("JetFlavourIdentifier",
        physicsDefinition = cms.bool(False),
        srcByReference = cms.InputTag("hltPatJetPartonAssociationLegacy")
    )
    process.hltPatJetPartonAssociationLegacy = cms.EDProducer("JetPartonMatcher",
        coneSizeToAssociate = cms.double(0.3),
        jets = cms.InputTag(pfjets),
        partons = cms.InputTag("hltPatJetPartonsLegacy")
    )
    process.hltPatJetPartonsLegacy = cms.EDProducer("PartonSelector",
        src = cms.InputTag("genParticles"),
        withLeptons = cms.bool(False)
    )
    process.hltPatJetGenJetMatch = cms.EDProducer("GenJetMatcher",
        checkCharge = cms.bool(False),
        matched = cms.InputTag("hltSlimmedGenJets"),
        maxDeltaR = cms.double(0.4),
        mcPdgId = cms.vint32(),
        mcStatus = cms.vint32(),
        resolveAmbiguities = cms.bool(True),
        resolveByMatchQuality = cms.bool(False),
        src = cms.InputTag(pfjets)
    )
    process.hltSlimmedGenJets = cms.EDProducer("PATGenJetSlimmer",
        clearDaughters = cms.bool(False),
        cut = cms.string('pt > 8'),
        cutLoose = cms.string(''),
        dropSpecific = cms.bool(False),
        nLoose = cms.uint32(0),
        packedGenParticles = cms.InputTag("hltPackedGenParticles"),
        src = cms.InputTag("ak4GenJetsNoNu")
    )
    process.hltPackedGenParticles = cms.EDProducer("PATPackedGenParticleProducer",
        inputCollection = cms.InputTag("hltPrunedGenParticlesWithStatusOne"),
        inputOriginal = cms.InputTag("genParticles"),
        map = cms.InputTag("hltPrunedGenParticles"),
        maxRapidity = cms.double(6)
    )
    process.hltPrunedGenParticlesWithStatusOne = cms.EDProducer("GenParticlePruner",
        select = cms.vstring(
            'drop  *',
            '++keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15',
            'drop   status == 2',
            'keep++ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)',
            'drop status == 1',
            'keep+ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)',
            'keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15',
            'keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16',
            '+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())',
            '+keep abs(pdgId) == 11 && status == 1 && (pt > 3 || isPromptFinalState())',
            'keep++ abs(pdgId) == 15',
            'drop  status > 30 && status < 70 ',
            'drop  pdgId == 21 && pt < 5',
            'drop   status == 2 && abs(pdgId) == 21',
            'keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 6 || abs(pdgId) == 37 ',
            'keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ',
            '+keep abs(pdgId) == 13 && status == 1',
            'keep (4 <= abs(pdgId) <= 5)',
            'keep (1 <= abs(pdgId) <= 3 || pdgId = 21) & (status = 2 || status = 11 || status = 71 || status = 72) && pt>5',
            'keep  abs(pdgId) == 323  && abs(eta) < 2.5 && pt > 1',
            'keep+ abs(pdgId) == 333',
            'keep+ abs(pdgId) == 9920443 || abs(pdgId) == 9042413 || abs(pdgId) == 9000443 || abs(pdgId) == 100541 || abs(pdgId) == 100543',
            'keep+ abs(pdgId) == 443 || abs(pdgId) == 100443 || abs(pdgId) == 10441 || abs(pdgId) == 20443 || abs(pdgId) == 445 || abs(pdgId) == 30443',
            'keep+ abs(pdgId) == 553 || abs(pdgId) == 100553 || abs(pdgId) == 200553 || abs(pdgId) == 10551 || abs(pdgId) == 20553 || abs(pdgId) == 555',
            'keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433',
            'keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543',
            'keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)',
            'keep (4900001 <= abs(pdgId) <= 4900991)',
            'keep (51 <= abs(pdgId) <= 53)',
            'keep pdgId = 2212',
            'keep status == 3 || ( 21 <= status <= 29) || ( 11 <= status <= 19)',
            'keep isHardProcess() || fromHardProcessFinalState() || fromHardProcessDecayed() || fromHardProcessBeforeFSR() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())',
            'keep    status == 1'
        ),
        src = cms.InputTag("genParticles")
    )
    process.hltPrunedGenParticles = cms.EDProducer("GenParticlePruner",
        select = cms.vstring(
            'drop  *',
            '++keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15',
            'drop   status == 2',
            'keep++ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)',
            'drop status == 1',
            'keep+ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)',
            'keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15',
            'keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16',
            '+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())',
            '+keep abs(pdgId) == 11 && status == 1 && (pt > 3 || isPromptFinalState())',
            'keep++ abs(pdgId) == 15',
            'drop  status > 30 && status < 70 ',
            'drop  pdgId == 21 && pt < 5',
            'drop   status == 2 && abs(pdgId) == 21',
            'keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 6 || abs(pdgId) == 37 ',
            'keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ',
            '+keep abs(pdgId) == 13 && status == 1',
            'keep (4 <= abs(pdgId) <= 5)',
            'keep (1 <= abs(pdgId) <= 3 || pdgId = 21) & (status = 2 || status = 11 || status = 71 || status = 72) && pt>5',
            'keep  abs(pdgId) == 323  && abs(eta) < 2.5 && pt > 1',
            'keep+ abs(pdgId) == 333',
            'keep+ abs(pdgId) == 9920443 || abs(pdgId) == 9042413 || abs(pdgId) == 9000443 || abs(pdgId) == 100541 || abs(pdgId) == 100543',
            'keep+ abs(pdgId) == 443 || abs(pdgId) == 100443 || abs(pdgId) == 10441 || abs(pdgId) == 20443 || abs(pdgId) == 445 || abs(pdgId) == 30443',
            'keep+ abs(pdgId) == 553 || abs(pdgId) == 100553 || abs(pdgId) == 200553 || abs(pdgId) == 10551 || abs(pdgId) == 20553 || abs(pdgId) == 555',
            'keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433',
            'keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543',
            'keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)',
            'keep (4900001 <= abs(pdgId) <= 4900991)',
            'keep (51 <= abs(pdgId) <= 53)',
            'keep pdgId = 2212',
            'keep status == 3 || ( 21 <= status <= 29) || ( 11 <= status <= 19)',
            'keep isHardProcess() || fromHardProcessFinalState() || fromHardProcessDecayed() || fromHardProcessBeforeFSR() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())'
        ),
        src = cms.InputTag("hltPrunedGenParticlesWithStatusOne")
    )
    process.hltPatJetPartonMatch = cms.EDProducer("MCMatcher",
        checkCharge = cms.bool(False),
        matched = cms.InputTag("hltPrunedGenParticles"),
        maxDPtRel = cms.double(3.0),
        maxDeltaR = cms.double(0.4),
        mcPdgId = cms.vint32(
            1, 2, 3, 4, 5,
            21
        ),
        mcStatus = cms.vint32(3, 23),
        resolveAmbiguities = cms.bool(True),
        resolveByMatchQuality = cms.bool(False),
        src = cms.InputTag(pfjets)
    )
    process.hltPatJetCharge = cms.EDProducer("JetChargeProducer",
        exp = cms.double(1.0),
        src = cms.InputTag("hltAk4JetTracksAssociatorAtVertexPF"),
        var = cms.string('Pt')
    )
    process.hltAk4JetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
        coneSize = cms.double(0.4),
        jets = cms.InputTag(pfjets),
        pvSrc = cms.InputTag(hltVertices),
        tracks = cms.InputTag(tracks),
        useAssigned = cms.bool(False)
    )
    process.hltPatJetCorrFactors = cms.EDProducer("JetCorrFactorsProducer",
        emf = cms.bool(False),
        extraJPTOffset = cms.string('L1FastJet'),
        flavorType = cms.string('J'),
        levels = cms.vstring(
            'L1FastJet',
            'L2Relative',
            'L3Absolute'
        ),
        payload = cms.string(payload),
        primaryVertices = cms.InputTag(hltVertices),
        rho = cms.InputTag(rho),
        src = cms.InputTag(pfjets),
        useNPV = cms.bool(True),
        useRho = cms.bool(True)
    )
    process.hltAk4JetID = cms.EDProducer("JetIDProducer",
        ebRecHitsColl = cms.InputTag(ecalRecHit,"EcalRecHitsEB"),
        eeRecHitsColl = cms.InputTag(ecalRecHit,"EcalRecHitsEE"),
        hbheRecHitsColl = cms.InputTag(hbhereco),
        hfRecHitsColl = cms.InputTag(hfreco),
        hoRecHitsColl = cms.InputTag(horeco),
        rpcRecHits = cms.InputTag(rpcRecHits),
        src = cms.InputTag(calojets),
        useRecHits = cms.bool(True)
    )
    process.hltPixelClusterTagInfos = cms.EDProducer("PixelClusterTagInfoProducer",
        addForward = cms.bool(True),
        hadronMass = cms.double(12.0),
        isPhase1 = cms.bool(True),
        jets = cms.InputTag(pfjets),
        maxJetEtaCut = cms.double(2.5),
        minAdcCount = cms.int32(-1),
        minJetPtCut = cms.double(100.0),
        pixelhit = cms.InputTag(siPixelClusters),
        vertices = cms.InputTag(hltVertices)
    )

    # process.hltPatJetCorrFactorsTransientCorrectedSlimmedDeepFlavour = cms.EDProducer("JetCorrFactorsProducer",
    #     emf = cms.bool(False),
    #     extraJPTOffset = cms.string('L1FastJet'),
    #     flavorType = cms.string('J'),
    #     levels = cms.vstring(
    #         'L1FastJet',
    #         'L2Relative',
    #         'L3Absolute'
    #     ),
    #     payload = cms.string(payload),
    #     primaryVertices = cms.InputTag(hltVerticesSlimmed),
    #     rho = cms.InputTag(rho),
    #     # src = cms.InputTag("hltUpdatedPatJetsSlimmedDeepFlavour"),
    #     src = cms.InputTag("hltPatJets"),
    #     useNPV = cms.bool(True),
    #     useRho = cms.bool(True)
    # )

    #### TAGGERS
    # run DeepFlavour for HLT
    process.hltPFDeepFlavourJetTags = cms.EDProducer("DeepFlavourONNXJetTagsProducer",
        flav_names = cms.vstring(
            'probb',
            'probbb',
            'problepb',
            'probc',
            'probuds',
            'probg'
        ),
        input_names = cms.vstring(
            'input_1',
            'input_2',
            'input_3',
            'input_4',
            'input_5'
        ),
        mightGet = cms.optional.untracked.vstring,
        model_path = cms.FileInPath('RecoBTag/Combined/data/DeepFlavourV03_10X_training/model.onnx'),
        output_names = cms.vstring('ID_pred/Softmax:0'),
        src = cms.InputTag("hltPFDeepFlavourTagInfos")
    )

    # process.hltPFDeepFlavourTagInfos = cms.EDProducer("DeepFlavourTagInfoProducer",
    #     candidates = cms.InputTag("hltPackedPFCandidates"),
    #     compute_probabilities = cms.bool(False),
    #     fallback_puppi_weight = cms.bool(False),
    #     fallback_vertex_association = cms.bool(False),
    #     flip = cms.bool(False),
    #     jet_radius = cms.double(0.4),
    #     # jets = cms.InputTag("hltUpdatedPatJetsSlimmedDeepFlavour"),
    #     # jets = cms.InputTag("hltPatJets"),
    #     jets = cms.InputTag(pfjets),
    #     max_jet_eta = cms.double(2.5),
    #     mightGet = cms.optional.untracked.vstring,
    #     min_candidate_pt = cms.double(0.95),
    #     min_jet_pt = cms.double(15),
    #     puppi_value_map = cms.InputTag(""),
    #     run_deepVertex = cms.bool(False),
    #     secondary_vertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
    #     shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagsInfos"),
    #     vertex_associator = cms.InputTag(""),
    #     vertices = cms.InputTag(hltVerticesSlimmed)
    # )

    process.hltPFDeepFlavourTagInfos = cms.EDProducer("DeepFlavourTagInfoProducer",
        # candidates = cms.InputTag("hltPackedPFCandidates"),
        candidates = cms.InputTag(particleFlow),
        compute_probabilities = cms.bool(False),
        fallback_puppi_weight = cms.bool(False),
        fallback_vertex_association = cms.bool(False),
        flip = cms.bool(False),
        jet_radius = cms.double(0.4),
        jets = cms.InputTag(pfjets),
        max_jet_eta = cms.double(2.5),
        mightGet = cms.optional.untracked.vstring,
        min_candidate_pt = cms.double(0.95),
        min_jet_pt = cms.double(15),
        puppi_value_map = cms.InputTag(puppi),
        run_deepVertex = cms.bool(False),
        secondary_vertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
        shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagsInfos"),
        vertex_associator = cms.InputTag("hltPrimaryVertexAssociation","original"),
        vertices = cms.InputTag(hltVertices)
    )


    process.candidateCombinedSecondaryVertexV2Computer = cms.ESProducer("CandidateCombinedSecondaryVertexESProducer",
        SoftLeptonFlip = cms.bool(False),
        calibrationRecords = cms.vstring(
            'CombinedSVIVFV2RecoVertex',
            'CombinedSVIVFV2PseudoVertex',
            'CombinedSVIVFV2NoVertex'
        ),
        categoryVariableName = cms.string('vertexCategory'),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        recordLabel = cms.string(''),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useCategories = cms.bool(True),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    )

    process.hltOfflineSlimmedPrimaryVertices = cms.EDProducer("PATVertexSlimmer",
        score = cms.InputTag("hltPrimaryVertexAssociation","original"),
        src = cms.InputTag(hltVertices)
    )

    # process.hltSlimmedSecondaryVertices = cms.EDProducer("PATSecondaryVertexSlimmer",
    #     packedPFCandidates = cms.InputTag("hltPackedPFCandidates"),
    #     src = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF")
    # )

    process.hltPrimaryVertexAssociation = cms.EDProducer("PFCandidatePrimaryVertexSorter",
        assignment = cms.PSet(
            maxDistanceToJetAxis = cms.double(0.07),
            maxDtSigForPrimaryAssignment = cms.double(3.0),
            maxDxyForJetAxisAssigment = cms.double(0.1),
            maxDxyForNotReconstructedPrimary = cms.double(0.01),
            maxDxySigForNotReconstructedPrimary = cms.double(2),
            maxDzErrorForPrimaryAssignment = cms.double(0.05),
            maxDzForJetAxisAssigment = cms.double(0.1),
            maxDzForPrimaryAssignment = cms.double(0.1),
            maxDzSigForPrimaryAssignment = cms.double(5.0),
            maxJetDeltaR = cms.double(0.5),
            minJetPt = cms.double(25),
            preferHighRanked = cms.bool(False),
            useTiming = cms.bool(False)
        ),
        jets = cms.InputTag(pfjets),
        particles = cms.InputTag(particleFlow),
        produceAssociationToOriginalVertices = cms.bool(True),
        produceNoPileUpCollection = cms.bool(False),
        producePileUpCollection = cms.bool(False),
        produceSortedVertices = cms.bool(False),
        qualityForPrimary = cms.int32(2),
        sorting = cms.PSet(

        ),
        usePVMET = cms.bool(True),
        vertices = cms.InputTag(hltVertices)
    )


    process.hltChargedHadronPFTrackIsolation = cms.EDProducer("ChargedHadronPFTrackIsolationProducer",
        mightGet = cms.optional.untracked.vstring,
        minRawCaloEnergy = cms.double(0.5),
        minTrackPt = cms.double(1),
        src = cms.InputTag(particleFlow)
    )


    # process.MC_NoFilter_PFBTagDeepCSV = cms.Path(process.HLTAK4PFJetsSequence+process.HLTBtagDeepCSVSequencePF)
    # process.MC_NoFilter_PFBTagDeepCSVPat = cms.Path(process.HLTAK4PFJetsSequence+process.HLTBtagDeepCSVSequencePFPat)

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
        *process.hltPixelClusterTagInfos
        *process.hltAk4JetTracksAssociatorAtVertexPF
        *process.hltPatJetCharge
        *process.hltPatJetCorrFactors

        *process.hltPrimaryVertexAssociation
        *process.hltOfflineSlimmedPrimaryVertices
        *process.hltChargedHadronPFTrackIsolation
        *process.hltPackedPFCandidates

        *process.hltPFDeepFlavourTagInfos
        *process.hltPFDeepFlavourJetTags

        *process.hltPatJets
        )

    if process.schedule_():
        # process.schedule.extend([process.MC_NoFilter_PFBTagDeepCSVPat])
        process.schedule.extend([process.MC_JetsMatchingPath])
        # process.schedule.extend([process.MC_NoFilter_PFBTagDeepCSV])

    return process
