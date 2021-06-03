import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak4PFClusterJets_cfi import ak4PFClusterJets as _ak4PFClusterJets
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJetsPuppi as _ak4PFJetsPuppi
from RecoJets.JetProducers.ak8PFJets_cfi import ak8PFJetsPuppi as _ak8PFJetsPuppi
from CommonTools.PileupAlgos.Puppi_cff import puppi as _puppi, puppiNoLep as _puppiNoLep
from JMETriggerAnalysis.Common.multiplicityValueProducerFromNestedCollectionEdmNewDetSetVectorSiPixelClusterDouble_cfi\
 import multiplicityValueProducerFromNestedCollectionEdmNewDetSetVectorSiPixelClusterDouble as _nSiPixelClusters


def addPaths_MC_JMEPFPuppiROI(process):

    process.hltPreMCJMEPFPuppiROI = cms.EDFilter('HLTPrescaler',
      L1GtReadoutRecordTag = cms.InputTag('hltGtStage2Digis'),
      offset = cms.uint32(0)
    )

    process.hltPixelClustersMultiplicity = _nSiPixelClusters.clone(src = 'hltSiPixelClusters', defaultValue = -1.)

    process.hltPFPuppiROI = _puppi.clone(
      candName = 'hltParticleFlowForBTag',
      vertexName = 'hltVerticesPFForBTag',
      usePUProxyValue = True,
      PUProxyValue = 'hltPixelClustersMultiplicity',
    )

    process.HLTPFPuppiSequenceROI = cms.Sequence(
        process.HLTDoCaloSequencePF
      + process.HLTL2muonrecoSequence
      + process.HLTL3muonrecoSequence
      + process.HLTTrackReconstructionForBTag
      + process.HLTParticleFlowSequenceForBTag
      + process.hltVerticesPFForBTag
      + process.hltPixelClustersMultiplicity
      + process.hltPFPuppiROI
    )

    ## AK4
    process.hltAK4PFPuppiJetsROI = _ak4PFJetsPuppi.clone(
      src = 'hltParticleFlowForBTag',
      srcWeights = 'hltPFPuppiROI',
      applyWeight = True,
    )

    process.hltAK4PFPuppiJetCorrectorL1ROI = cms.EDProducer('L1FastjetCorrectorProducer',
      algorithm = cms.string('AK4PFPuppiHLT'),
      level = cms.string('L1FastJet'),
      srcRho = cms.InputTag('hltFixedGridRhoFastjetAllForBTag'),
    )

    process.hltAK4PFPuppiJetCorrectorL2ROI = cms.EDProducer('LXXXCorrectorProducer',
      algorithm = cms.string('AK4PFPuppiHLT'),
      level = cms.string('L2Relative')
    )

    process.hltAK4PFPuppiJetCorrectorL3ROI = cms.EDProducer('LXXXCorrectorProducer',
      algorithm = cms.string('AK4PFPuppiHLT'),
      level = cms.string('L3Absolute')
    )

    process.hltAK4PFPuppiJetCorrectorROI = cms.EDProducer('ChainedJetCorrectorProducer',
      correctors = cms.VInputTag(
        'hltAK4PFPuppiJetCorrectorL1ROI',
        'hltAK4PFPuppiJetCorrectorL2ROI',
        'hltAK4PFPuppiJetCorrectorL3ROI',
      ),
    )

    process.hltAK4PFPuppiJetsCorrectedROI = cms.EDProducer('CorrectedPFJetProducer',
                                                           src = cms.InputTag('hltAK4PFPuppiJetsROI'),
                                                           correctors = cms.VInputTag('hltAK4PFPuppiJetCorrectorROI'),
    )

    process.HLTAK4PFPuppiJetsSequenceROI = cms.Sequence(
        process.hltAK4PFPuppiJetsROI
      + process.hltAK4PFPuppiJetCorrectorL1ROI
      + process.hltAK4PFPuppiJetCorrectorL2ROI
      + process.hltAK4PFPuppiJetCorrectorL3ROI
      + process.hltAK4PFPuppiJetCorrectorROI
      + process.hltAK4PFPuppiJetsCorrectedROI
    )


    ## Modifications to PUPPI parameters
    for mod_i in [process.hltPFPuppiROI]:
      for algo_idx in range(len(mod_i.algos)):
        if len(mod_i.algos[algo_idx].MinNeutralPt) != len(mod_i.algos[algo_idx].MinNeutralPtSlope):
          raise RuntimeError('instance of PuppiProducer is misconfigured:\n\n'+str(mod_i)+' = '+mod_i.dumpPython())

        for algoReg_idx in range(len(mod_i.algos[algo_idx].MinNeutralPt)):
          mod_i.algos[algo_idx].MinNeutralPt[algoReg_idx] += 2.56 * mod_i.algos[algo_idx].MinNeutralPtSlope[algoReg_idx]
          mod_i.algos[algo_idx].MinNeutralPtSlope[algoReg_idx] *= 0.00271

    ## Path
    process.MC_JMEPFPuppiROI_v1 = cms.Path(
        process.HLTBeginSequence
        + process.hltPreMCJMEPFPuppiROI
        ## AK4 Jets
        + process.HLTAK4PFJetsSequenceForBTag # New
        + process.HLTPFPuppiSequenceROI
        + process.HLTAK4PFPuppiJetsSequenceROI
    )

    if process.schedule_():
       process.schedule_().append(process.MC_JMEPFPuppiROI_v1)

    return process


def addPaths_PFJetsForBtag(process):

    process.hltPreMCROIPFBTagDeepCSV = cms.EDFilter("HLTPrescaler",
                                                    L1GtReadoutRecordTag = cms.InputTag("hltGtStage2Digis"),
                                                    offset = cms.uint32(0)
    )

    process.hltLightPFTracksForBTag = process.hltLightPFTracks.clone(
        TkColList = cms.VInputTag("hltMergedTracksForBTag"),
    )


    process.hltParticleFlowBlockForBTag = process.hltParticleFlowBlock.clone(
        elementImporters = cms.VPSet(
            cms.PSet(
                DPtOverPtCuts_byTrackAlgo = cms.vdouble(
                    0.5, 0.5, 0.5, 0.5, 0.5,
                    0.5
                ),
                NHitCuts_byTrackAlgo = cms.vuint32(
                    3, 3, 3, 3, 3,
                    3
                ),
                cleanBadConvertedBrems = cms.bool(False),
                importerName = cms.string('GeneralTracksImporter'),
                muonMaxDPtOPt = cms.double(1.0),
                muonSrc = cms.InputTag("hltMuons"),
                source = cms.InputTag("hltLightPFTracksForBTag"),
                trackQuality = cms.string('highPurity'),
                useIterativeTracking = cms.bool(False)
            ),
        )
    )

    process.hltParticleFlowForBTag = process.hltParticleFlow.clone(
        blocks = cms.InputTag("hltParticleFlowBlockForBTag"),
        vertexCollection = cms.InputTag("hltPixelVertices")
    )


    process.HLTParticleFlowSequenceForBTag = cms.Sequence(
        process.HLTPreshowerSequence
        + process.hltParticleFlowRecHitECALUnseeded
        + process.hltParticleFlowRecHitHBHE
        + process.hltParticleFlowRecHitHF
        + process.hltParticleFlowRecHitPSUnseeded
        + process.hltParticleFlowClusterECALUncorrectedUnseeded
        + process.hltParticleFlowClusterPSUnseeded
        + process.hltParticleFlowClusterECALUnseeded
        + process.hltParticleFlowClusterHBHE
        + process.hltParticleFlowClusterHCAL
        + process.hltParticleFlowClusterHF
        + process.hltLightPFTracksForBTag
        + process.hltParticleFlowBlockForBTag
        + process.hltParticleFlowForBTag)

    process.hltAK4PFJetsForBTag = process.hltAK4PFJets.clone(
        src = cms.InputTag("hltParticleFlowForBTag"),
        srcPVs = cms.InputTag("hltPixelVertices"),
    )

    process.hltAK4PFJetsLooseIDForBTag = process.hltAK4PFJetsLooseID.clone(
        jetsInput = cms.InputTag("hltAK4PFJetsForBTag"),
    )


    process.hltAK4PFJetsTightIDForBTag = process.hltAK4PFJetsTightID.clone(
        jetsInput = cms.InputTag("hltAK4PFJetsForBTag"),
    )

    process.HLTAK4PFJetsReconstructionSequenceForBTag = cms.Sequence(
        process.HLTL2muonrecoSequence
        + process.HLTL3muonrecoSequence
        + process.HLTTrackReconstructionForBTag
        + process.HLTParticleFlowSequenceForBTag 
        + process.hltAK4PFJetsForBTag
        + process.hltAK4PFJetsLooseIDForBTag
        + process.hltAK4PFJetsTightIDForBTag
    )

    process.hltFixedGridRhoFastjetAllForBTag = process.hltFixedGridRhoFastjetAll.clone(
        pfCandidatesTag = cms.InputTag("hltParticleFlowForBTag")
    )


    process.hltAK4PFFastJetCorrectorForBTag = process.hltAK4PFFastJetCorrector.clone(
        srcRho = cms.InputTag("hltFixedGridRhoFastjetAllForBTag")
    )

    process.hltAK4PFCorrectorForBTag = cms.EDProducer("ChainedJetCorrectorProducer",
                                                      correctors = cms.VInputTag("hltAK4PFFastJetCorrectorForBTag", "hltAK4PFRelativeCorrector", "hltAK4PFAbsoluteCorrector", "hltAK4PFResidualCorrector")
    )


    process.HLTAK4PFCorrectorProducersSequenceForBTag = cms.Sequence(
        process.hltAK4PFFastJetCorrectorForBTag
        + process.hltAK4PFRelativeCorrector
        + process.hltAK4PFAbsoluteCorrector
        + process.hltAK4PFResidualCorrector
        + process.hltAK4PFCorrectorForBTag
    )


    process.hltAK4PFJetsCorrectedForBTag = cms.EDProducer("CorrectedPFJetProducer",
                                                          correctors = cms.VInputTag("hltAK4PFCorrectorForBTag"),
                                                          src = cms.InputTag("hltAK4PFJetsForBTag")
    )


    process.hltAK4PFJetsLooseIDCorrectedForBTag = cms.EDProducer("CorrectedPFJetProducer",
                                                                 correctors = cms.VInputTag("hltAK4PFCorrectorForBTag"),
                                                                 src = cms.InputTag("hltAK4PFJetsLooseIDForBTag")
    )


    process.hltAK4PFJetsTightIDCorrectedForBTag = cms.EDProducer("CorrectedPFJetProducer",
                                                                 correctors = cms.VInputTag("hltAK4PFCorrectorForBTag"),
                                                                 src = cms.InputTag("hltAK4PFJetsTightIDForBTag")
    )

    process.HLTAK4PFJetsCorrectionSequenceForBTag = cms.Sequence(
        process.hltFixedGridRhoFastjetAllForBTag
        + process.HLTAK4PFCorrectorProducersSequenceForBTag
        + process.hltAK4PFJetsCorrectedForBTag
        + process.hltAK4PFJetsLooseIDCorrectedForBTag
        + process.hltAK4PFJetsTightIDCorrectedForBTag
    )

    process.HLTAK4PFJetsSequenceForBTag = cms.Sequence(
        process.HLTPreAK4PFJetsRecoSequence
        + process.HLTAK4PFJetsReconstructionSequenceForBTag
        + process.HLTAK4PFJetsCorrectionSequenceForBTag
    )

    process.hltVerticesPFForBTag = process.hltVerticesPF.clone(
        TrackLabel = cms.InputTag("hltMergedTracksForBTag"),
    )

    process.hltVerticesPFSelectorForBTag = process.hltVerticesPFSelector.clone(
        filterParams = cms.PSet(
            maxRho = cms.double(2.0),
            maxZ = cms.double(24.0),
            minNdof = cms.double(4.0),
            pvSrc = cms.InputTag("hltVerticesPFForBTag")
        ),
        src = cms.InputTag("hltVerticesPFForBTag")
    )

    process.hltVerticesPFFilterForBTag = process.hltVerticesPFFilter.clone(
        src = cms.InputTag("hltVerticesPFSelectorForBTag")
    )

    process.hltPFJetForBtagSelectorForBTag = process.hltPFJetForBtagSelector.clone(
        # inputTag = cms.InputTag("hltAK4PFJetsCorrected"),
        inputTag = cms.InputTag("hltAK4PFJetsCorrectedForBTag"), 
    )

    process.hltPFJetForBtagROI = process.hltPFJetForBtag.clone(
        HLTObject = cms.InputTag("hltPFJetForBtagSelectorForBTag"),
    )

    process.hltDeepBLifetimeTagInfosPFROI = process.hltDeepBLifetimeTagInfosPF.clone(
        candidates = cms.InputTag("hltParticleFlowForBTag"),
        jets = cms.InputTag("hltPFJetForBtagROI"),
        primaryVertex = cms.InputTag("hltVerticesPFFilterForBTag"),
    )


    process.hltDeepInclusiveVertexFinderPFROI = process.hltDeepInclusiveVertexFinderPF.clone(
        primaryVertices = cms.InputTag("hltVerticesPFFilterForBTag"),
        tracks = cms.InputTag("hltParticleFlowForBTag"),
    )

    process.hltDeepInclusiveSecondaryVerticesPFROI = process.hltDeepInclusiveSecondaryVerticesPF.clone(
        secondaryVertices = cms.InputTag("hltDeepInclusiveVertexFinderPFROI")
    )

    process.hltDeepTrackVertexArbitratorPFROI = process.hltDeepTrackVertexArbitratorPF.clone(
        primaryVertices = cms.InputTag("hltVerticesPFFilterForBTag"),
        secondaryVertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPFROI"),
        tracks = cms.InputTag("hltParticleFlowForBTag")
    )
    
    process.hltDeepInclusiveMergedVerticesPFROI = process.hltDeepInclusiveMergedVerticesPF.clone(
        secondaryVertices = cms.InputTag("hltDeepTrackVertexArbitratorPFROI")
    )

    process.hltDeepSecondaryVertexTagInfosPFROI = process.hltDeepSecondaryVertexTagInfosPF.clone(
        extSVCollection = cms.InputTag("hltDeepInclusiveMergedVerticesPFROI"),
        trackIPTagInfos = cms.InputTag("hltDeepBLifetimeTagInfosPFROI"),
    )

    process.hltDeepCombinedSecondaryVertexBJetTagsInfosROI = process.hltDeepCombinedSecondaryVertexBJetTagsInfos.clone(
        svTagInfos = cms.InputTag("hltDeepSecondaryVertexTagInfosPFROI")
    )

    process.hltDeepCombinedSecondaryVertexBJetTagsPFROI = process.hltDeepCombinedSecondaryVertexBJetTagsPF.clone(
        src = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsInfosROI"),
    )

    process.HLTBtagDeepCSVSequencePFROI = cms.Sequence(
        process.hltVerticesPFForBTag
        + process.hltVerticesPFSelectorForBTag
        + process.hltVerticesPFFilterForBTag
        + process.hltPFJetForBtagSelectorForBTag
        + process.hltPFJetForBtagROI
        + process.hltDeepBLifetimeTagInfosPFROI
        + process.hltDeepInclusiveVertexFinderPFROI
        + process.hltDeepInclusiveSecondaryVerticesPFROI
        + process.hltDeepTrackVertexArbitratorPFROI
        + process.hltDeepInclusiveMergedVerticesPFROI
        + process.hltDeepSecondaryVertexTagInfosPFROI
        + process.hltDeepCombinedSecondaryVertexBJetTagsInfosROI
        + process.hltDeepCombinedSecondaryVertexBJetTagsPFROI
    )

    process.hltBTagPFDeepCSV4p06SingleROI = process.hltBTagPFDeepCSV4p06Single.clone(
        JetTags = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsPFROI","probb"),
        Jets = cms.InputTag("hltPFJetForBtag"),
    )


    process.MC_ROIPFBTagDeepCSV_v10 = cms.Path(
        process.HLTBeginSequence
        + process.hltPreMCROIPFBTagDeepCSV
        + process.HLTAK4PFJetsSequenceForBTag
        + process.HLTBtagDeepCSVSequencePFROI
        + process.hltBTagPFDeepCSV4p06SingleROI
        + process.HLTEndSequence)
        

    if process.schedule_():
       process.schedule_().append(process.MC_ROIPFBTagDeepCSV_v10)

    return process
