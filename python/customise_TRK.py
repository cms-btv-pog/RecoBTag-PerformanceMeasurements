import FWCore.ParameterSet.Config as cms

def customisePFForPixelTracks(process, tracksToUse = "hltPixelTracks"):
     process.hltPFMuonMerging.TrackProducers = cms.VInputTag("hltIterL3MuonTracks", tracksToUse)
     process.hltPFMuonMerging.selectedTrackQuals = cms.VInputTag("hltIterL3MuonTracks", tracksToUse)
     process.hltParticleFlowBlock.elementImporters = cms.VPSet(
             cms.PSet(
                 # DPtOverPtCuts_byTrackAlgo = cms.vdouble(
                 #     0.5, 0.5, 0.5, 0.5, 0.5,
                 #     0.5
                 # ),
                 DPtOverPtCuts_byTrackAlgo = cms.vdouble(
                     5.0, 5.0, 5.0, 5.0, 5.0,
                     5.0
                 ),
                 NHitCuts_byTrackAlgo = cms.vuint32(
                     3, 3, 3, 3, 3,
                     3
                 ),
                 cleanBadConvertedBrems = cms.bool(False),
                 importerName = cms.string('GeneralTracksImporter'),
                 muonMaxDPtOPt = cms.double(1.0),
                 muonSrc = cms.InputTag("hltMuons"),
                 source = cms.InputTag("hltLightPFTracks"),
                 # trackQuality = cms.string('highPurity'),
                 trackQuality = cms.string('any'),
                 useIterativeTracking = cms.bool(False)
             ),
             cms.PSet(
                 BCtoPFCMap = cms.InputTag(""),
                 importerName = cms.string('ECALClusterImporter'),
                 source = cms.InputTag("hltParticleFlowClusterECALUnseeded")
             ),
             cms.PSet(
                 importerName = cms.string('GenericClusterImporter'),
                 source = cms.InputTag("hltParticleFlowClusterHCAL")
             ),
             cms.PSet(
                 importerName = cms.string('GenericClusterImporter'),
                 source = cms.InputTag("hltParticleFlowClusterHF")
             )
         )
     process.hltMuonLinks = process.hltPixelOnlyMuonLinks.clone(
         InclusiveTrackerTrackCollection = cms.InputTag("hltPFMuonMerging")
     )
     process.hltMuons = process.hltPixelOnlyMuons.clone(
         TrackExtractorPSet = cms.PSet(
             BeamSpotLabel = cms.InputTag("hltOnlineBeamSpot"),
             BeamlineOption = cms.string('BeamSpotFromEvent'),
             Chi2Ndof_Max = cms.double(1e+64),
             Chi2Prob_Min = cms.double(-1.0),
             ComponentName = cms.string('TrackExtractor'),
             DR_Max = cms.double(1.0),
             DR_Veto = cms.double(0.01),
             DepositLabel = cms.untracked.string(''),
             Diff_r = cms.double(0.1),
             Diff_z = cms.double(0.2),
             NHits_Min = cms.uint32(0),
             Pt_Min = cms.double(-1.0),
             inputTrackCollection = cms.InputTag("hltPFMuonMerging")
         ),
         inputCollectionLabels = cms.VInputTag("hltPFMuonMerging", "hltMuonLinks", "hltL2Muons")
     )
     process.HLTTrackReconstructionForPF = cms.Sequence(
         process.HLTDoLocalPixelSequence+
         process.HLTRecopixelvertexingSequence+
         process.hltPFMuonMerging+
         process.hltMuonLinks+
         process.hltMuons)

     return process

def customisePFForPixelTracksCleaned(process, tracksToUse = "hltPixelTracksCleanForBTag", vertex="hltTrimmedPixelVertices", nVertices = 4):
    process.hltPixelTracksCleanForBTag = cms.EDProducer("TrackWithVertexSelector",
        copyExtras = cms.untracked.bool(False),
        copyTrajectories = cms.untracked.bool(False),
        d0Max = cms.double(999.0),
        dzMax = cms.double(999.0),
        etaMax = cms.double(5.0),
        etaMin = cms.double(0.0),
        nSigmaDtVertex = cms.double(0.0),
        nVertices = cms.uint32(nVertices),
        normalizedChi2 = cms.double(999999.0),
        numberOfLostHits = cms.uint32(999),
        numberOfValidHits = cms.uint32(0),
        numberOfValidHitsForGood = cms.uint32(999),
        numberOfValidPixelHits = cms.uint32(3),
        numberOfValidPixelHitsForGood = cms.uint32(999),
        ptErrorCut = cms.double(5.0),
        ptMax = cms.double(500.0),
        # ptMin = cms.double(0.3),
        ptMin = cms.double(0.8),
        quality = cms.string('any'),
        rhoVtx = cms.double(0.2),
        rhoVtxScale = cms.double(1.0),
        rhoVtxSig = cms.double(999.0),
        src = cms.InputTag("hltPixelTracks"),
        timeResosTag = cms.InputTag(""),
        timesTag = cms.InputTag(""),
        useVtx = cms.bool(True),
        vertexTag = cms.InputTag(vertex),
        vtxFallback = cms.bool(True),
        zetaVtx = cms.double(0.3),
        zetaVtxScale = cms.double(1.0),
        zetaVtxSig = cms.double(999.0)
    )
    process.hltPFMuonMerging.TrackProducers = cms.VInputTag("hltIterL3MuonTracks", tracksToUse)
    process.hltPFMuonMerging.selectedTrackQuals = cms.VInputTag("hltIterL3MuonTracks", tracksToUse)
    process.hltParticleFlowBlock.elementImporters = cms.VPSet(
        cms.PSet(
             # DPtOverPtCuts_byTrackAlgo = cms.vdouble(
             #     0.5, 0.5, 0.5, 0.5, 0.5,
             #     0.5
             # ),
             DPtOverPtCuts_byTrackAlgo = cms.vdouble(
                 5.0, 5.0, 5.0, 5.0, 5.0,
                 5.0
             ),
             NHitCuts_byTrackAlgo = cms.vuint32(
                 3, 3, 3, 3, 3,
                 3
             ),
             # NHitCuts_byTrackAlgo = cms.vuint32(
             #     0, 0, 0, 0, 0,
             #     0
             # ),
             cleanBadConvertedBrems = cms.bool(False),
             importerName = cms.string('GeneralTracksImporter'),
             muonMaxDPtOPt = cms.double(1.0),
             muonSrc = cms.InputTag("hltMuons"),
             source = cms.InputTag("hltLightPFTracks"),
             # trackQuality = cms.string('highPurity'),
             trackQuality = cms.string('any'),
             useIterativeTracking = cms.bool(False)
        ),
        cms.PSet(
            BCtoPFCMap = cms.InputTag(""),
            importerName = cms.string('ECALClusterImporter'),
            source = cms.InputTag("hltParticleFlowClusterECALUnseeded")
        ),
        cms.PSet(
            importerName = cms.string('GenericClusterImporter'),
            source = cms.InputTag("hltParticleFlowClusterHCAL")
        ),
        cms.PSet(
            importerName = cms.string('GenericClusterImporter'),
            source = cms.InputTag("hltParticleFlowClusterHF")
        )
    )
    process.hltMuonLinks = process.hltPixelOnlyMuonLinks.clone(
        InclusiveTrackerTrackCollection = cms.InputTag("hltPFMuonMerging")
    )
    process.hltMuons = process.hltPixelOnlyMuons.clone(
        TrackExtractorPSet = cms.PSet(
            BeamSpotLabel = cms.InputTag("hltOnlineBeamSpot"),
            BeamlineOption = cms.string('BeamSpotFromEvent'),
            Chi2Ndof_Max = cms.double(1e+64),
            Chi2Prob_Min = cms.double(-1.0),
            ComponentName = cms.string('TrackExtractor'),
            DR_Max = cms.double(1.0),
            DR_Veto = cms.double(0.01),
            DepositLabel = cms.untracked.string(''),
            Diff_r = cms.double(0.1),
            Diff_z = cms.double(0.2),
            NHits_Min = cms.uint32(0),
            Pt_Min = cms.double(-1.0),
            inputTrackCollection = cms.InputTag("hltPFMuonMerging")
        ),
        inputCollectionLabels = cms.VInputTag("hltPFMuonMerging", "hltMuonLinks", "hltL2Muons")
    )
    process.HLTTrackReconstructionForPF = cms.Sequence(
        process.HLTDoLocalPixelSequence+
        process.HLTRecopixelvertexingSequence+
        process.hltPixelTracksCleanForBTag+
        process.hltPFMuonMerging+
        process.hltMuonLinks+
        process.hltMuons
    )
    process.hltVerticesPF.TkFilterParameters.minSiliconLayersWithHits = cms.int32(0)
    return process

def customizeVertices(process):
    process.hltPixelVerticesPFSelector = cms.EDFilter("PrimaryVertexObjectFilter",
        filterParams = cms.PSet(
            maxRho = cms.double(2.0),
            maxZ = cms.double(24.0),
            minNdof = cms.double(4.0),
            pvSrc = cms.InputTag("hltPixelVertices")
        ),
        src = cms.InputTag("hltPixelVertices")
    )

    process.hltPixelVerticesPFFilter = cms.EDFilter("VertexSelector",
        cut = cms.string('!isFake'),
        filter = cms.bool(True),
        src = cms.InputTag("hltPixelVerticesPFSelector")
    )

    process.HLTBtagDeepCSVSequencePFPat = cms.Sequence(
        process.hltPixelVerticesPFSelector
        + process.hltPixelVerticesPFFilter
        + process.hltVerticesPF
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
    process.HLTBtagDeepCSVSequencePFPuppiPat = cms.Sequence(
        process.hltPixelVerticesPFSelector
        + process.hltPixelVerticesPFFilter
        + process.hltVerticesPF
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

    process.hltAk4JetTracksAssociatorAtVertexPF.pvSrc = cms.InputTag("hltPixelVerticesPFFilter")
    process.hltAk4JetTracksAssociatorAtVertexPFPuppi.pvSrc = cms.InputTag("hltPixelVerticesPFFilter")

    process.hltDeepBLifetimePFPatTagInfos.primaryVertex = cms.InputTag("hltPixelVerticesPFFilter")
    process.hltDeepBLifetimePFPuppiPatTagInfos.primaryVertex = cms.InputTag("hltPixelVerticesPFFilter")
    process.hltDeepBLifetimeTagInfosPF.primaryVertex = cms.InputTag("hltPixelVerticesPFFilter")
    process.hltDeepInclusiveVertexFinderPF.primaryVertices = cms.InputTag("hltPixelVerticesPFFilter")
    process.hltDeepTrackVertexArbitratorPF.primaryVertices = cms.InputTag("hltPixelVerticesPFFilter")
    process.hltPFDeepFlavourTagInfos.vertices = cms.InputTag("hltPixelVerticesPFFilter")
    process.hltPFPuppiDeepFlavourTagInfos.vertices = cms.InputTag("hltPixelVerticesPFFilter")

    return process

def customizeVertices2(process):
    process.hltVerticesPFFilter.src = cms.InputTag("hltPixelVertices")
    process.hltPFPuppi.vertexName = cms.InputTag("hltPixelVertices")
    process.hltPFPuppiNoLep.vertexName = cms.InputTag("hltPixelVertices")

    return process

def customizeMinHitsAndPt(process, doForPat=False):
    # ptValue = 0.8
    ptValue = 0.9

    process.hltPixelTracksCleanForBTag.ptMin = cms.double(ptValue)

    process.hltDeepInclusiveVertexFinderPF.minHits = cms.uint32(0)
    process.hltDeepBLifetimeTagInfosPF.minimumNumberOfHits = cms.int32(0)
    process.hltImpactParameterTagInfos.minimumNumberOfHits = cms.int32(0)
    process.hltDeepTrackVertexArbitratorPF.trackMinLayers = cms.int32(0)
    process.hltDeepSecondaryVertexTagInfosPF.trackSelection.totalHitsMin = cms.uint32(2)
    process.hltDeepCombinedSecondaryVertexBJetTagsInfos.computer.trackPseudoSelection.totalHitsMin = cms.uint32(0)
    process.hltDeepCombinedSecondaryVertexBJetTagsInfos.computer.trackSelection.totalHitsMin = cms.uint32(2)
    process.hltVerticesPF.TkFilterParameters.minPt = cms.double(ptValue)
    process.hltDeepInclusiveVertexFinderPF.minPt = cms.double(ptValue)
    process.hltDeepTrackVertexArbitratorPF.trackMinPt = cms.double(ptValue)
    process.hltDeepCombinedSecondaryVertexBJetTagsInfos.computer.trackSelection.ptMin = cms.double(ptValue)
    if doForPat:
        process.hltDeepBLifetimePFPatTagInfos.minimumNumberOfHits = cms.int32(0)
        process.hltImpactParameterPatTagInfos.minimumNumberOfHits = cms.int32(0)
        process.hltDeepSecondaryVertexPFPatTagInfos.trackSelection.totalHitsMin = cms.uint32(2)
        process.hltDeepCombinedSecondaryVertexBJetPatTagInfos.computer.trackPseudoSelection.totalHitsMin = cms.uint32(0)
        process.hltDeepCombinedSecondaryVertexBJetPatTagInfos.computer.trackSelection.totalHitsMin = cms.uint32(2)
        process.hltDeepCombinedSecondaryVertexBJetPatTagInfos.computer.trackSelection.ptMin = cms.double(ptValue)

    return process

def addDeepJet(process, doPF=True, doPuppi=False):
    from PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi import primaryVertexAssociation
    from RecoBTag.FeatureTools.pfDeepFlavourTagInfos_cfi import pfDeepFlavourTagInfos
    from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags

    if doPF:
        process.hltPrimaryVertexAssociation = primaryVertexAssociation.clone(
            jets = cms.InputTag("hltPFJetForBtag"),
            particles = cms.InputTag("hltParticleFlow" ),
            vertices = cms.InputTag("hltVerticesPFFilter"),
        )
        process.hltPFDeepFlavourTagInfos = pfDeepFlavourTagInfos.clone(
            candidates = cms.InputTag("hltParticleFlow" ),
            jets = cms.InputTag("hltPFJetForBtag"),
            fallback_puppi_weight = cms.bool(True),
            puppi_value_map = cms.InputTag(""),
            secondary_vertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
            # shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBJetPatTagInfos"),
            shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsInfos"),
            vertex_associator = cms.InputTag("hltPrimaryVertexAssociation","original"),
            vertices = cms.InputTag("hltVerticesPFFilter")
        )
        process.hltPFDeepFlavourJetTags = pfDeepFlavourJetTags.clone(
            src = cms.InputTag("hltPFDeepFlavourTagInfos")
        )

        process.HLTBtagDeepJetSequencePF = cms.Sequence(
            process.hltVerticesPF
            +process.hltVerticesPFSelector
            +process.hltVerticesPFFilter
            +process.hltPFJetForBtagSelector
            +process.hltPFJetForBtag
            +process.hltDeepBLifetimeTagInfosPF
            +process.hltDeepInclusiveVertexFinderPF
            +process.hltDeepInclusiveSecondaryVerticesPF
            +process.hltDeepTrackVertexArbitratorPF
            +process.hltDeepInclusiveMergedVerticesPF
            +process.hltDeepSecondaryVertexTagInfosPF
            +process.hltDeepCombinedSecondaryVertexBJetTagsInfos
            +process.hltPrimaryVertexAssociation
            +process.hltPFDeepFlavourTagInfos
            +process.hltPFDeepFlavourJetTags
        )

        process.hltPreMCPFBTagDeepJet = process.hltPreMCPFBTagDeepCSV.clone()

        process.hltBTagPFDeepJet4p06Single = process.hltBTagPFDeepCSV4p06Single.clone(
            JetTags = cms.InputTag("hltPFDeepFlavourJetTags","probb"),
            Jets = cms.InputTag("hltPFJetForBtag"),
            MaxTag = cms.double(999999.0),
            MinJets = cms.int32(1),
            MinTag = cms.double(0.25),
            TriggerType = cms.int32(86),
            saveTags = cms.bool(True)
        )

        process.MC_PFBTagDeepJet = cms.Path(
            process.HLTBeginSequence
            +process.hltPreMCPFBTagDeepCSV
            +process.HLTAK4PFJetsSequence
            +process.HLTBtagDeepJetSequencePF
            +process.hltBTagPFDeepJet4p06Single
            +process.HLTEndSequence
        )

    if doPuppi:
        process.hltPFPuppiDeepFlavourJetTags = pfDeepFlavourJetTags.clone(
            src = cms.InputTag("hltPFPuppiDeepFlavourTagInfos")
        )
        process.hltPFPuppiDeepFlavourTagInfos = pfDeepFlavourTagInfos.clone(
            candidates = cms.InputTag("hltParticleFlow" ),
            jets = cms.InputTag("hltAK4PFPuppiJets"),
            puppi_value_map = cms.InputTag("hltPFPuppi"),
            secondary_vertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
            shallow_tag_infos = cms.InputTag("hltDeepCombinedSecondaryVertexBPuppiJetPatTagInfos"),
            vertex_associator = cms.InputTag("hltPrimaryVertexAssociationPuppi","original"),
            vertices = cms.InputTag("hltVerticesPFFilter")
        )
        process.hltPrimaryVertexAssociationPuppi = primaryVertexAssociation.clone(
            jets = cms.InputTag("hltAK4PFPuppiJets"),
            particles = cms.InputTag("hltParticleFlow" ),
            vertices = cms.InputTag("hltVerticesPFFilter"),
        )
    if process.schedule_():
        if doPF: process.schedule.extend([process.MC_PFBTagDeepJet])
        # if doPuppi: process.schedule.extend([process.MC_PuppiJetsMatchingPath])
    return process
