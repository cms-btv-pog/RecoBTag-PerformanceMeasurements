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
      #+ process.HLTL2muonrecoSequence
      #+ process.HLTL3muonrecoSequence
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
        + process.HLTPFPuppiSequenceROI
        + process.HLTAK4PFPuppiJetsSequenceROI
    )

    if process.schedule_():
       process.schedule_().append(process.MC_JMEPFPuppiROI_v1)

    return process
