#ifndef JETINFOBRANCHES_H
#define JETINFOBRANCHES_H

#include <TTree.h>
#include "RecoBTag/PerformanceMeasurements/interface/VariableParser.h"

const UInt_t nMaxJets_ = 10000;
const UInt_t nMaxTrk_  = 100000;
const UInt_t nMaxMuons_= 10000;
const UInt_t nMaxElectrons_= 10000;
const UInt_t nMaxSVs_= 10000;
const UInt_t nMaxLeptons_=10000;

class JetInfoBranches {

  public :

    int   nJet;
    float Jet_pt[nMaxJets_];
    float Jet_uncorrpt[nMaxJets_];
    float Jet_genpt[nMaxJets_];
    float Jet_residual[nMaxJets_];
    float Jet_area[nMaxJets_];
    float Jet_jes[nMaxJets_];
    float Jet_eta[nMaxJets_];
    float Jet_phi[nMaxJets_];
    float Jet_mass[nMaxJets_];
    float Jet_Ip1N[nMaxJets_];
    float Jet_Ip1P[nMaxJets_];
    float Jet_Ip2N[nMaxJets_];
    float Jet_Ip2P[nMaxJets_];
    float Jet_Ip3N[nMaxJets_];
    float Jet_Ip3P[nMaxJets_];
    float Jet_Ip4N[nMaxJets_];
    float Jet_Ip4P[nMaxJets_];
    float Jet_Mass4N[nMaxJets_];
    float Jet_Mass4P[nMaxJets_];

    float Jet_DeepFlavourBDisc[nMaxJets_]   ;
    float Jet_DeepFlavourCvsLDisc[nMaxJets_];
    float Jet_DeepFlavourCvsBDisc[nMaxJets_];
    float Jet_DeepFlavourB[nMaxJets_];
    float Jet_DeepFlavourBB[nMaxJets_];
    float Jet_DeepFlavourLEPB[nMaxJets_];
    float Jet_DeepFlavourC[nMaxJets_];
    float Jet_DeepFlavourUDS[nMaxJets_];
    float Jet_DeepFlavourG[nMaxJets_];

    float Jet_DeepFlavourBDiscN[nMaxJets_]   ;
    float Jet_DeepFlavourCvsLDiscN[nMaxJets_];
    float Jet_DeepFlavourCvsBDiscN[nMaxJets_];

    float Jet_DeepFlavourPrunedBDisc[nMaxJets_]   ;
    float Jet_DeepFlavourPrunedCvsLDisc[nMaxJets_];
    float Jet_DeepFlavourPrunedCvsBDisc[nMaxJets_];
    float Jet_DeepFlavourPrunedB[nMaxJets_];
    float Jet_DeepFlavourPrunedBB[nMaxJets_];
    float Jet_DeepFlavourPrunedLEPB[nMaxJets_];
    float Jet_DeepFlavourPrunedC[nMaxJets_];
    float Jet_DeepFlavourPrunedUDS[nMaxJets_];
    float Jet_DeepFlavourPrunedG[nMaxJets_];

    float Jet_DeepFlavourPrunedBDiscN[nMaxJets_]   ;
    float Jet_DeepFlavourPrunedCvsLDiscN[nMaxJets_];
    float Jet_DeepFlavourPrunedCvsBDiscN[nMaxJets_];


    float Jet_DeepCSVBDisc[nMaxJets_]  ;
    float Jet_DeepCSVBDiscN[nMaxJets_]  ;
    float Jet_DeepCSVBDiscP[nMaxJets_]  ;
    float Jet_DeepCSVCvsLDisc[nMaxJets_]  ;
    float Jet_DeepCSVCvsLDiscN[nMaxJets_]  ;
    float Jet_DeepCSVCvsLDiscP[nMaxJets_]  ;
    float Jet_DeepCSVCvsBDisc[nMaxJets_]  ;
    float Jet_DeepCSVCvsBDiscN[nMaxJets_]  ;
    float Jet_DeepCSVCvsBDiscP[nMaxJets_]  ;
    float Jet_DeepCSVb[nMaxJets_]  ;
    float Jet_DeepCSVc[nMaxJets_]  ;
    float Jet_DeepCSVl[nMaxJets_]  ;
    float Jet_DeepCSVbb[nMaxJets_] ;
    float Jet_DeepCSVcc[nMaxJets_] ;
    float Jet_DeepCSVbN[nMaxJets_] ;
    float Jet_DeepCSVcN[nMaxJets_] ;
    float Jet_DeepCSVlN[nMaxJets_] ;
    float Jet_DeepCSVbbN[nMaxJets_];
    float Jet_DeepCSVccN[nMaxJets_];
    float Jet_DeepCSVbP[nMaxJets_] ;
    float Jet_DeepCSVcP[nMaxJets_] ;
    float Jet_DeepCSVlP[nMaxJets_] ;
    float Jet_DeepCSVbbP[nMaxJets_];
    float Jet_DeepCSVccP[nMaxJets_];

    float Jet_ProbaN[nMaxJets_];
    float Jet_ProbaP[nMaxJets_];
    float Jet_Proba[nMaxJets_];
    float Jet_BprobN[nMaxJets_];
    float Jet_BprobP[nMaxJets_];
    float Jet_Bprob[nMaxJets_];
    float Jet_SvxN[nMaxJets_];
    float Jet_Svx[nMaxJets_];
    float Jet_SvxNHP[nMaxJets_];
    float Jet_SvxHP[nMaxJets_];
    float Jet_CombSvxN[nMaxJets_];
    float Jet_CombSvxP[nMaxJets_];
    float Jet_CombSvx[nMaxJets_];
    float Jet_CombIVF[nMaxJets_];
    float Jet_CombIVF_P[nMaxJets_];
    float Jet_CombIVF_N[nMaxJets_];
    float Jet_SoftMuN[nMaxJets_];
    float Jet_SoftMuP[nMaxJets_];
    float Jet_SoftMu[nMaxJets_];
    float Jet_SoftElN[nMaxJets_];
    float Jet_SoftElP[nMaxJets_];
    float Jet_SoftEl[nMaxJets_];
    float Jet_DoubleSV[nMaxJets_];
    float Jet_cMVA[nMaxJets_];
    float Jet_cMVAv2[nMaxJets_];
    float Jet_cMVAv2N[nMaxJets_];
    float Jet_cMVAv2P[nMaxJets_];
    //New test variables for AK4 jets: to be cleaned up in the future
    float Jet_trackSip2dSig_AboveBottom_0[nMaxJets_];
    float Jet_trackSip2dSig_AboveBottom_1[nMaxJets_];
    //End new test variables
    int   Jet_hist1[nMaxJets_];
    int   Jet_hist2[nMaxJets_];
    int   Jet_hist3[nMaxJets_];
    int   Jet_histJet[nMaxJets_];
    int   Jet_histSvx[nMaxJets_];
    int   Jet_ntracks[nMaxJets_];
    int   Jet_nseltracks[nMaxJets_];
    int   Jet_flavour[nMaxJets_];	
    int   Jet_flavourCleaned[nMaxJets_];
    int   Jet_partonFlavour[nMaxJets_];
    int   Jet_hadronFlavour[nMaxJets_];
    int   Jet_partonid[nMaxJets_];
    int   Jet_nbHadrons[nMaxJets_];
    int   Jet_ncHadrons[nMaxJets_];
    int   Jet_nFirstTrack[nMaxJets_];
    int   Jet_nLastTrack[nMaxJets_];
    int   Jet_nFirstTrackTruth[nMaxJets_];
    int   Jet_nLastTrackTruth[nMaxJets_];
    int   Jet_nFirstSV[nMaxJets_];
    int   Jet_nLastSV[nMaxJets_];
    int   Jet_SV_multi[nMaxJets_];
    int   Jet_nFirstTrkInc[nMaxJets_];
    int   Jet_nLastTrkInc[nMaxJets_];
    int   Jet_nSM[nMaxJets_];
    int   Jet_nFirstSM[nMaxJets_];
    int   Jet_nLastSM[nMaxJets_];
    int   Jet_nSE[nMaxJets_];
    int   Jet_nFirstSE[nMaxJets_];
    int   Jet_nLastSE[nMaxJets_];
    int   Jet_looseID[nMaxJets_];
    int   Jet_tightID[nMaxJets_];
    int   Jet_FatJetIdx[nMaxJets_];
    float Jet_ptSoftDrop[nMaxJets_];
    float Jet_etaSoftDrop[nMaxJets_];
    float Jet_phiSoftDrop[nMaxJets_];
    float Jet_massSoftDrop[nMaxJets_];
    float Jet_jecF0SoftDrop[nMaxJets_];
    float Jet_ptPruned[nMaxJets_];
    float Jet_etaPruned[nMaxJets_];
    float Jet_phiPruned[nMaxJets_];
    float Jet_massPruned[nMaxJets_];
    float Jet_jecF0Pruned[nMaxJets_];
    float Jet_tau1[nMaxJets_];
    float Jet_tau2[nMaxJets_];
    float Jet_tauAxis1_px[nMaxJets_];
    float Jet_tauAxis1_py[nMaxJets_];
    float Jet_tauAxis1_pz[nMaxJets_];
    float Jet_tauAxis2_px[nMaxJets_];
    float Jet_tauAxis2_py[nMaxJets_];
    float Jet_tauAxis2_pz[nMaxJets_];
    float Jet_z_ratio[nMaxJets_];
    float Jet_nTracks_fat[nMaxJets_];
    float Jet_nSV_fat[nMaxJets_];
    float Jet_trackSip3dSig_3[nMaxJets_];
    float Jet_trackSip3dSig_2[nMaxJets_];
    float Jet_trackSip3dSig_1[nMaxJets_];
    float Jet_trackSip3dSig_0[nMaxJets_];
    float Jet_trackSip2dSigAboveCharm_0[nMaxJets_];
    float Jet_trackSip2dSigAboveCharm_1[nMaxJets_];
    float Jet_trackSip2dSigAboveBottom_0[nMaxJets_];
    float Jet_trackSip2dSigAboveBottom_1[nMaxJets_];
    float Jet_tau1_trackSip3dSig_0[nMaxJets_];
    float Jet_tau1_trackSip3dSig_1[nMaxJets_];
    float Jet_tau2_trackSip3dSig_0[nMaxJets_];
    float Jet_tau2_trackSip3dSig_1[nMaxJets_];
    float Jet_tau1_trackEtaRel_0[nMaxJets_];
    float Jet_tau1_trackEtaRel_1[nMaxJets_];
    float Jet_tau1_trackEtaRel_2[nMaxJets_];
    float Jet_tau2_trackEtaRel_0[nMaxJets_];
    float Jet_tau2_trackEtaRel_1[nMaxJets_];
    float Jet_tau2_trackEtaRel_2[nMaxJets_];
    float Jet_tau1_nSecondaryVertices[nMaxJets_];
    float Jet_tau2_nSecondaryVertices[nMaxJets_];
    float Jet_tau1_flightDistance2dSig[nMaxJets_];
    float Jet_tau2_flightDistance2dSig[nMaxJets_];
    float Jet_tau1_vertexDeltaR[nMaxJets_];
    float Jet_tau2_vertexDeltaR[nMaxJets_];
    float Jet_tau1_vertexEnergyRatio[nMaxJets_];
    float Jet_tau2_vertexEnergyRatio[nMaxJets_];
    float Jet_tau1_vertexMass[nMaxJets_];
    float Jet_tau2_vertexMass[nMaxJets_];
    float Jet_tau1_vertexMass_corrected[nMaxJets_];
    float Jet_tau2_vertexMass_corrected[nMaxJets_];
    float Jet_tau1_vertexNTracks[nMaxJets_];
    float Jet_tau2_vertexNTracks[nMaxJets_];    
    int   Jet_nFirstTrkTagVar[nMaxJets_];
    int   Jet_nLastTrkTagVar[nMaxJets_];
    int   Jet_nFirstSVTagVar[nMaxJets_];
    int   Jet_nLastSVTagVar[nMaxJets_];
    int   Jet_nFirstTrkTagVarCSV[nMaxJets_];
    int   Jet_nLastTrkTagVarCSV[nMaxJets_];
    int   Jet_nFirstTrkEtaRelTagVarCSV[nMaxJets_];
    int   Jet_nLastTrkEtaRelTagVarCSV[nMaxJets_];
    float Jet_DeepDoubleB[nMaxJets_];

    int   nTrack;
    float Track_dxy[nMaxTrk_];
    float Track_dz[nMaxTrk_];
    float Track_dxyError[nMaxTrk_];
    float Track_dzError[nMaxTrk_];
    int   Track_sign2D[nMaxTrk_];
    int   Track_sign3D[nMaxTrk_];
    float Track_LongIP[nMaxTrk_];
    float Track_length[nMaxTrk_];
    float Track_dist[nMaxTrk_];
    float Track_IP2D[nMaxTrk_];
    float Track_IP2Dsig[nMaxTrk_];
    float Track_IP2Derr[nMaxTrk_];
    float Track_IP[nMaxTrk_];
    float Track_IPsig[nMaxTrk_];
    float Track_IPerr[nMaxTrk_];
    float Track_Proba[nMaxTrk_];
    float Track_p[nMaxTrk_];
    float Track_pt[nMaxTrk_];
    float Track_eta[nMaxTrk_];
    float Track_phi[nMaxTrk_];
    float Track_chi2[nMaxTrk_];
    int   Track_charge[nMaxTrk_];
    int   Track_history[nMaxTrk_];
    int   Track_nHitStrip[nMaxTrk_];
    int   Track_nHitPixel[nMaxTrk_];
    int   Track_nHitAll[nMaxTrk_];
    int   Track_nHitTIB[nMaxTrk_];
    int   Track_nHitTID[nMaxTrk_];
    int   Track_nHitTOB[nMaxTrk_];
    int   Track_nHitTEC[nMaxTrk_];
    int   Track_nHitPXB[nMaxTrk_];
    int   Track_nHitPXF[nMaxTrk_];
    int   Track_isHitL1[nMaxTrk_];
    int   Track_PV[nMaxTrk_];
    int   Track_SV[nMaxTrk_];
    int   Track_isfromSV[nMaxTrk_];
    float Track_PVweight[nMaxTrk_];
    float Track_SVweight[nMaxTrk_];
    int   Track_isfromV0[nMaxTrk_];
    float Track_lengthTau[nMaxTrk_];
    float Track_distTau[nMaxTrk_];    
    int   Track_category[nMaxTrk_];
    
    float Track_TPAssociationQuality[nMaxTrk_];
    int	  Track_idxMatchedTP[nMaxTrk_];
    int   nTrackTruth;
    int   TrackTruth_idxMatchedTrack[nMaxTrk_];
    float TrackTruth_p[nMaxTrk_];
    float TrackTruth_pt[nMaxTrk_];
    float TrackTruth_eta[nMaxTrk_];
    float TrackTruth_phi[nMaxTrk_];
    int   TrackTruth_charge[nMaxTrk_];
    int   TrackTruth_pdgid[nMaxTrk_];
    float TrackTruth_dxy[nMaxTrk_];
    float TrackTruth_dz[nMaxTrk_];
    int   TrackTruth_nHitAll[nMaxTrk_];
    int   TrackTruth_nHitPixel[nMaxTrk_];
    int   TrackTruth_nHitStrip[nMaxTrk_];

    int   nTrkInc;
    float TrkInc_pt[nMaxTrk_];
    float TrkInc_eta[nMaxTrk_];
    float TrkInc_phi[nMaxTrk_];
    float TrkInc_ptrel[nMaxTrk_];
    float TrkInc_IPsig[nMaxTrk_];
    float TrkInc_IP[nMaxTrk_];

    int   nPFElectron;
    int   PFElectron_IdxJet[nMaxElectrons_];
    float PFElectron_pt[nMaxElectrons_];
    float PFElectron_eta[nMaxElectrons_];
    float PFElectron_phi[nMaxElectrons_];
    float PFElectron_ptrel[nMaxElectrons_];
    float PFElectron_ratio[nMaxElectrons_];
    float PFElectron_ratioRel[nMaxElectrons_];
    float PFElectron_deltaR[nMaxElectrons_];
    float PFElectron_IP[nMaxElectrons_];
    float PFElectron_IP2D[nMaxElectrons_];
    float PFElectron_mva_e_pi[nMaxElectrons_];

    int   nPFMuon;
    int   PFMuon_IdxJet[nMaxMuons_];
    int   PFMuon_nMuHit[nMaxMuons_];
    int   PFMuon_nTkHit[nMaxMuons_];
    int   PFMuon_nPixHit[nMaxMuons_];
    int   PFMuon_nOutHit[nMaxMuons_];
    int   PFMuon_nTkLwM[nMaxMuons_];
    int   PFMuon_nPixLwM[nMaxMuons_];
    int   PFMuon_nMatched[nMaxMuons_];
    float PFMuon_chi2[nMaxMuons_];
    float PFMuon_chi2Tk[nMaxMuons_];
    int   PFMuon_isGlobal[nMaxMuons_];
    int   PFMuon_hist[nMaxMuons_];
    float PFMuon_pt[nMaxMuons_];
    float PFMuon_eta[nMaxMuons_];
    float PFMuon_phi[nMaxMuons_];
    float PFMuon_ptrel[nMaxMuons_];
    float PFMuon_ratio[nMaxMuons_];
    float PFMuon_ratioRel[nMaxMuons_];
    float PFMuon_deltaR[nMaxMuons_];
    float PFMuon_IP[nMaxMuons_];
    float PFMuon_IP2D[nMaxMuons_];
    float PFMuon_IPsig[nMaxMuons_];
    float PFMuon_IP2Dsig[nMaxMuons_];
    float PFMuon_dz[nMaxMuons_];
    int   PFMuon_GoodQuality[nMaxMuons_];

    int   nSV;
    float SV_x[nMaxSVs_];
    float SV_y[nMaxSVs_];
    float SV_z[nMaxSVs_];
    float SV_ex[nMaxSVs_];
    float SV_ey[nMaxSVs_];
    float SV_ez[nMaxSVs_];
    float SV_chi2[nMaxSVs_];
    float SV_ndf[nMaxSVs_];
    float SV_flight[nMaxSVs_];
    float SV_flightErr[nMaxSVs_];
    float SV_deltaR_jet[nMaxSVs_];
    float SV_deltaR_sum_jet[nMaxSVs_];
    float SV_deltaR_sum_dir[nMaxSVs_];
    float SV_vtx_pt[nMaxSVs_];
    float SV_flight2D[nMaxSVs_];
    float SV_flight2DErr[nMaxSVs_];
    float SV_totCharge[nMaxSVs_];
    float SV_vtxDistJetAxis[nMaxSVs_];
    int   SV_nTrk[nMaxSVs_];
    float SV_mass[nMaxSVs_];
    float SV_vtx_eta[nMaxSVs_];
    float SV_vtx_phi[nMaxSVs_];
    float SV_EnergyRatio[nMaxSVs_];
    float SV_dir_x[nMaxSVs_];
    float SV_dir_y[nMaxSVs_];
    float SV_dir_z[nMaxSVs_];

    // TagInfo TaggingVariables
    // per jet
    // run 2 id: https://twiki.cern.ch/CMS/JetID13TeVRun2017
    int TagVar_jetNTracks[nMaxJets_];
    int TagVar_jetNSecondaryVertices[nMaxJets_];
    float TagVar_chargedHadronEnergyFraction[nMaxJets_];
    float TagVar_neutralHadronEnergyFraction[nMaxJets_];
    float TagVar_photonEnergyFraction[nMaxJets_];
    float TagVar_electronEnergyFraction[nMaxJets_];
    float TagVar_muonEnergyFraction[nMaxJets_];
    int TagVar_chargedHadronMultiplicity[nMaxJets_];
    int TagVar_neutralHadronMultiplicity[nMaxJets_];
    int TagVar_photonMultiplicity[nMaxJets_];
    int TagVar_electronMultiplicity[nMaxJets_];
    int TagVar_muonMultiplicity[nMaxJets_];
    float TagVar_neutralEmEnergyFraction[nMaxJets_];
    int TagVar_numberOfDaughters[nMaxJets_];
    int TagVar_chargedMultiplicity[nMaxJets_];
    float TagVar_chargedEmEnergyFraction[nMaxJets_];
    int TagVar_neutralMultiplicity[nMaxJets_];

    // per jet per track
    int   nTrkTagVar;
    float TagVar_trackMomentum[nMaxTrk_];
    float TagVar_trackEta[nMaxTrk_];
    float TagVar_trackPhi[nMaxTrk_];
    float TagVar_trackPtRel[nMaxTrk_];
    float TagVar_trackPPar[nMaxTrk_];
    float TagVar_trackEtaRel[nMaxTrk_];
    float TagVar_trackDeltaR[nMaxTrk_];
    float TagVar_trackPtRatio[nMaxTrk_];
    float TagVar_trackPParRatio[nMaxTrk_];
    float TagVar_trackSip2dVal[nMaxTrk_];
    float TagVar_trackSip2dSig[nMaxTrk_];
    float TagVar_trackSip3dVal[nMaxTrk_];
    float TagVar_trackSip3dSig[nMaxTrk_];
    float TagVar_trackDecayLenVal[nMaxTrk_];
    float TagVar_trackDecayLenSig[nMaxTrk_];
    float TagVar_trackJetDistVal[nMaxTrk_];
    float TagVar_trackJetDistSig[nMaxTrk_];
    float TagVar_trackChi2[nMaxTrk_];
    float TagVar_trackNTotalHits[nMaxTrk_];
    float TagVar_trackNPixelHits[nMaxTrk_];
    // per jet per secondary vertex
    int   nSVTagVar;
    float TagVar_vertexMass[nMaxSVs_];
    float TagVar_vertexNTracks[nMaxSVs_];
    float TagVar_vertexJetDeltaR[nMaxSVs_];
    float TagVar_flightDistance2dVal[nMaxSVs_];
    float TagVar_flightDistance2dSig[nMaxSVs_];
    float TagVar_flightDistance3dVal[nMaxSVs_];
    float TagVar_flightDistance3dSig[nMaxSVs_];

    //DeepFlavour input features selection
    float DeepFlavourInput_charged_Sip3dVal[nMaxJets_];
    float DeepFlavourInput_charged_Sip3dSig[nMaxJets_];
    float DeepFlavourInput_charged_quality[ nMaxJets_];
    float DeepFlavourInput_charged_chi2[    nMaxJets_];

    float DeepFlavourInput_neutral_drminsv[nMaxJets_];
    float DeepFlavourInput_neutral_hadFrac[nMaxJets_];
    float DeepFlavourInput_neutral_ptrel[  nMaxJets_];

    float DeepFlavourInput_sv_d3d[     nMaxJets_];
    float DeepFlavourInput_sv_d3dsig[  nMaxJets_];
    float DeepFlavourInput_sv_normchi2[nMaxJets_];

    // CSV TaggingVariables
    // per jet
    float TagVarCSV_trackJetPt[nMaxJets_];
    float TagVarCSV_jetNTracks[nMaxJets_];
    float TagVarCSV_jetNTracksEtaRel[nMaxJets_];
    float TagVarCSV_trackSumJetEtRatio[nMaxJets_];
    float TagVarCSV_trackSumJetDeltaR[nMaxJets_];
    float TagVarCSV_trackSip2dValAboveCharm[nMaxJets_];
    float TagVarCSV_trackSip2dSigAboveCharm[nMaxJets_];
    float TagVarCSV_trackSip3dValAboveCharm[nMaxJets_];
    float TagVarCSV_trackSip3dSigAboveCharm[nMaxJets_];
    float TagVarCSV_vertexCategory[nMaxJets_];
    float TagVarCSV_jetNSecondaryVertices[nMaxJets_];
    float TagVarCSV_vertexMass[nMaxJets_];
    float TagVarCSV_vertexNTracks[nMaxJets_];
    float TagVarCSV_vertexEnergyRatio[nMaxJets_];
    float TagVarCSV_vertexJetDeltaR[nMaxJets_];
    float TagVarCSV_flightDistance2dVal[nMaxJets_];
    float TagVarCSV_flightDistance2dSig[nMaxJets_];
    float TagVarCSV_flightDistance3dVal[nMaxJets_];
    float TagVarCSV_flightDistance3dSig[nMaxJets_];
    // per jet per track
    int   nTrkTagVarCSV;
    int   nTrkEtaRelTagVarCSV;
    float TagVarCSV_trackMomentum[nMaxTrk_];
    float TagVarCSV_trackEta[nMaxTrk_];
    float TagVarCSV_trackPhi[nMaxTrk_];
    float TagVarCSV_trackPtRel[nMaxTrk_];
    float TagVarCSV_trackPPar[nMaxTrk_];
    float TagVarCSV_trackDeltaR[nMaxTrk_];
    float TagVarCSV_trackPtRatio[nMaxTrk_];
    float TagVarCSV_trackPParRatio[nMaxTrk_];
    float TagVarCSV_trackSip2dVal[nMaxTrk_];
    float TagVarCSV_trackSip2dSig[nMaxTrk_];
    float TagVarCSV_trackSip3dVal[nMaxTrk_];
    float TagVarCSV_trackSip3dSig[nMaxTrk_];
    float TagVarCSV_trackDecayLenVal[nMaxTrk_];
    float TagVarCSV_trackDecayLenSig[nMaxTrk_];
    float TagVarCSV_trackJetDistVal[nMaxTrk_];
    float TagVarCSV_trackJetDistSig[nMaxTrk_];
    float TagVarCSV_trackEtaRel[nMaxTrk_];

    //C tag TaggingVariables
    int   Jet_nFirstTrkCTagVar[nMaxJets_];
    int   Jet_nLastTrkCTagVar[nMaxJets_];
    int   Jet_nFirstTrkEtaRelCTagVar[nMaxJets_];
    int   Jet_nLastTrkEtaRelCTagVar[nMaxJets_];
    int   Jet_nFirstLepCTagVar[nMaxJets_];
    int   Jet_nLastLepCTagVar[nMaxJets_];
    //per jet    
    float CTag_Jet_CvsB[nMaxJets_];
    float CTag_Jet_CvsBN[nMaxJets_];
    float CTag_Jet_CvsBP[nMaxJets_];
    float CTag_Jet_CvsL[nMaxJets_];
    float CTag_Jet_CvsLN[nMaxJets_];
    float CTag_Jet_CvsLP[nMaxJets_];
    float CTag_jetNTracks[nMaxJets_];
    float CTag_jetNTracksEtaRel[nMaxJets_];
    float CTag_jetNLeptons[nMaxJets_]; 
    float CTag_trackSumJetEtRatio[nMaxJets_];
    float CTag_trackSumJetDeltaR[nMaxJets_];
    float CTag_trackSip2dSigAboveCharm[nMaxJets_];
    float CTag_trackSip3dSigAboveCharm[nMaxJets_];
    float CTag_vertexCategory[nMaxJets_];
    float CTag_jetNSecondaryVertices[nMaxJets_];
    float CTag_vertexMass[nMaxJets_];
    float CTag_vertexNTracks[nMaxJets_];
    float CTag_vertexEnergyRatio[nMaxJets_];
    float CTag_vertexJetDeltaR[nMaxJets_];
    float CTag_flightDistance2dSig[nMaxJets_];
    float CTag_flightDistance3dSig[nMaxJets_];
    float CTag_massVertexEnergyFraction[nMaxJets_];
    float CTag_vertexBoostOverSqrtJetPt[nMaxJets_];
    float CTag_vertexLeptonCategory[nMaxJets_];
    //per jet per track
    int   nTrkCTagVar;
    int   nTrkEtaRelCTagVar;
    float CTag_trackPtRel[nMaxTrk_];
    float CTag_trackPPar[nMaxTrk_];
    float CTag_trackDeltaR[nMaxTrk_];
    float CTag_trackPtRatio[nMaxTrk_];
    float CTag_trackPParRatio[nMaxTrk_];
    float CTag_trackSip2dSig[nMaxTrk_];
    float CTag_trackSip3dSig[nMaxTrk_];
    float CTag_trackDecayLenVal[nMaxTrk_];
    float CTag_trackJetDistVal[nMaxTrk_];
    float CTag_trackEtaRel[nMaxTrk_];
    //per jet per lepton
    int   nLeptons;
    float CTag_leptonPtRel[nMaxLeptons_];
    float CTag_leptonSip3d[nMaxLeptons_];
    float CTag_leptonDeltaR[nMaxLeptons_];
    float CTag_leptonRatioRel[nMaxLeptons_];
    float CTag_leptonEtaRel[nMaxLeptons_];
    float CTag_leptonRatio[nMaxLeptons_];

    void RegisterBranches(TTree *tree, VariableParser variableParser, std::string name="") {
      if(name!="") name += ".";
      if(variableParser.isToBeStored(name+"nJet")) tree->Branch((name+"nJet").c_str(),            &nJet           ,(name+"nJet/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_pt")) tree->Branch((name+"Jet_pt").c_str(),          Jet_pt          ,(name+"Jet_pt["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_uncorrpt")) tree->Branch((name+"Jet_uncorrpt").c_str(),    Jet_uncorrpt    ,(name+"Jet_uncorrpt["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_genpt")) tree->Branch((name+"Jet_genpt").c_str(),       Jet_genpt       ,(name+"Jet_genpt["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_residual")) tree->Branch((name+"Jet_residual").c_str(),    Jet_residual    ,(name+"Jet_residual["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_area")) tree->Branch((name+"Jet_area").c_str(),        Jet_area        ,(name+"Jet_area["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_jes")) tree->Branch((name+"Jet_jes").c_str(),         Jet_jes         ,(name+"Jet_jes["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_eta")) tree->Branch((name+"Jet_eta").c_str(),         Jet_eta         ,(name+"Jet_eta["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_phi")) tree->Branch((name+"Jet_phi").c_str(),         Jet_phi         ,(name+"Jet_phi["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_mass")) tree->Branch((name+"Jet_mass").c_str(),        Jet_mass        ,(name+"Jet_mass["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_ntracks")) tree->Branch((name+"Jet_ntracks").c_str(),     Jet_ntracks     ,(name+"Jet_ntracks["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nseltracks")) tree->Branch((name+"Jet_nseltracks").c_str(),  Jet_nseltracks  ,(name+"Jet_nseltracks["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_flavour")) tree->Branch((name+"Jet_flavour").c_str(),     Jet_flavour     ,(name+"Jet_flavour["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_flavourCleaned")) tree->Branch((name+"Jet_flavourCleaned").c_str(), Jet_flavourCleaned, (name+"Jet_flavourCleaned["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_partonFlavour")) tree->Branch((name+"Jet_partonFlavour").c_str(), Jet_partonFlavour, (name+"Jet_partonFlavour["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_hadronFlavour")) tree->Branch((name+"Jet_hadronFlavour").c_str(), Jet_hadronFlavour, (name+"Jet_hadronFlavour["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nbHadrons")) tree->Branch((name+"Jet_nbHadrons").c_str(),   Jet_nbHadrons   ,(name+"Jet_nbHadrons["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_ncHadrons")) tree->Branch((name+"Jet_ncHadrons").c_str(),   Jet_ncHadrons   ,(name+"Jet_ncHadrons["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_Ip2N")) tree->Branch((name+"Jet_Ip2N").c_str(),        Jet_Ip2N        ,(name+"Jet_Ip2N["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_Ip2P")) tree->Branch((name+"Jet_Ip2P").c_str(),        Jet_Ip2P        ,(name+"Jet_Ip2P["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_Ip3N")) tree->Branch((name+"Jet_Ip3N").c_str(),        Jet_Ip3N        ,(name+"Jet_Ip3N["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_Ip3P")) tree->Branch((name+"Jet_Ip3P").c_str(),        Jet_Ip3P        ,(name+"Jet_Ip3P["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_ProbaN")) tree->Branch((name+"Jet_ProbaN").c_str(),      Jet_ProbaN     ,(name+"Jet_ProbaN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_ProbaP")) tree->Branch((name+"Jet_ProbaP").c_str(),      Jet_ProbaP     ,(name+"Jet_ProbaP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_Proba")) tree->Branch((name+"Jet_Proba").c_str(),       Jet_Proba      ,(name+"Jet_Proba["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_BprobN")) tree->Branch((name+"Jet_BprobN").c_str(),      Jet_BprobN     ,(name+"Jet_BprobN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_BprobP")) tree->Branch((name+"Jet_BprobP").c_str(),      Jet_BprobP     ,(name+"Jet_BprobP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_Bprob")) tree->Branch((name+"Jet_Bprob").c_str(),       Jet_Bprob      ,(name+"Jet_Bprob["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_SvxN")) tree->Branch((name+"Jet_SvxN").c_str(),        Jet_SvxN       ,(name+"Jet_SvxN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_Svx")) tree->Branch((name+"Jet_Svx").c_str(),         Jet_Svx        ,(name+"Jet_Svx["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_SvxNHP")) tree->Branch((name+"Jet_SvxNHP").c_str(),      Jet_SvxNHP     ,(name+"Jet_SvxNHP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_SvxHP")) tree->Branch((name+"Jet_SvxHP").c_str(),       Jet_SvxHP      ,(name+"Jet_SvxHP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_CombSvxN")) tree->Branch((name+"Jet_CombSvxN").c_str(),    Jet_CombSvxN   ,(name+"Jet_CombSvxN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_CombSvxP")) tree->Branch((name+"Jet_CombSvxP").c_str(),    Jet_CombSvxP   ,(name+"Jet_CombSvxP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_CombSvx")) tree->Branch((name+"Jet_CombSvx").c_str(),     Jet_CombSvx    ,(name+"Jet_CombSvx["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_CombIVF")) tree->Branch((name+"Jet_CombIVF").c_str(),     Jet_CombIVF     ,(name+"Jet_CombIVF["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_CombIVF_P")) tree->Branch((name+"Jet_CombIVF_P").c_str(),   Jet_CombIVF_P   ,(name+"Jet_CombIVF_P["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_CombIVF_N")) tree->Branch((name+"Jet_CombIVF_N").c_str(),   Jet_CombIVF_N   ,(name+"Jet_CombIVF_N["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_SoftMuN")) tree->Branch((name+"Jet_SoftMuN").c_str(),     Jet_SoftMuN     ,(name+"Jet_SoftMuN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_SoftMuP")) tree->Branch((name+"Jet_SoftMuP").c_str(),     Jet_SoftMuP     ,(name+"Jet_SoftMuP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_SoftMu")) tree->Branch((name+"Jet_SoftMu").c_str(),      Jet_SoftMu      ,(name+"Jet_SoftMu["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_SoftElN")) tree->Branch((name+"Jet_SoftElN").c_str(),     Jet_SoftElN     ,(name+"Jet_SoftElN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_SoftElP")) tree->Branch((name+"Jet_SoftElP").c_str(),     Jet_SoftElP     ,(name+"Jet_SoftElP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_SoftEl")) tree->Branch((name+"Jet_SoftEl").c_str(),      Jet_SoftEl      ,(name+"Jet_SoftEl["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_cMVA")) tree->Branch((name+"Jet_cMVA").c_str(),    Jet_cMVA    ,(name+"Jet_cMVA["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_cMVAv2")) tree->Branch((name+"Jet_cMVAv2").c_str(),    Jet_cMVAv2    ,(name+"Jet_cMVAv2["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_cMVAv2N")) tree->Branch((name+"Jet_cMVAv2N").c_str(),    Jet_cMVAv2N    ,(name+"Jet_cMVAv2N["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_cMVAv2P")) tree->Branch((name+"Jet_cMVAv2P").c_str(),    Jet_cMVAv2P    ,(name+"Jet_cMVAv2P["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_hist1")) tree->Branch((name+"Jet_hist1").c_str(),       Jet_hist1       ,(name+"Jet_hist1["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_hist2")) tree->Branch((name+"Jet_hist2").c_str(),       Jet_hist2       ,(name+"Jet_hist2["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_hist3")) tree->Branch((name+"Jet_hist3").c_str(),       Jet_hist3       ,(name+"Jet_hist3["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_histJet")) tree->Branch((name+"Jet_histJet").c_str(),     Jet_histJet     ,(name+"Jet_histJet["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_histSvx")) tree->Branch((name+"Jet_histSvx").c_str(),     Jet_histSvx     ,(name+"Jet_histSvx["+name+"nJet]/I").c_str());

      if(variableParser.isToBeStored(name+"Jet_SV_multi")) tree->Branch((name+"Jet_SV_multi").c_str(),    Jet_SV_multi      ,(name+"Jet_SV_multi["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nSM")) tree->Branch((name+"Jet_nSM").c_str(),         Jet_nSM         ,(name+"Jet_nSM["+name+"nJet]/I").c_str()      );
      if(variableParser.isToBeStored(name+"Jet_nSE")) tree->Branch((name+"Jet_nSE").c_str(),         Jet_nSE         ,(name+"Jet_nSE["+name+"nJet]/I").c_str()      );

      if(variableParser.isToBeStored(name+"Jet_looseID")) tree->Branch((name+"Jet_looseID").c_str(),      Jet_looseID  ,(name+"Jet_looseID["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_tightID")) tree->Branch((name+"Jet_tightID").c_str(),      Jet_tightID  ,(name+"Jet_tightID["+name+"nJet]/I").c_str());

      //test variables for AK4 jets: clean up in the future
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSig_AboveBottom_0")) tree->Branch((name+"Jet_trackSip2dSig_AboveBottom_0").c_str(),     Jet_trackSip2dSig_AboveBottom_0     ,(name+"Jet_trackSip2dSig_AboveBottom_0["+name+"nJet]/F").c_str()            );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSig_AboveBottom_1")) tree->Branch((name+"Jet_trackSip2dSig_AboveBottom_1").c_str(),     Jet_trackSip2dSig_AboveBottom_1     ,(name+"Jet_trackSip2dSig_AboveBottom_1["+name+"nJet]/F").c_str()            );

      if(variableParser.isToBeStored(name+"Jet_DeepFlavourBDisc")) tree->Branch((name+"Jet_DeepFlavourBDisc"   ).c_str(),  Jet_DeepFlavourBDisc    , (name+"Jet_DeepFlavourBDisc["+name+"nJet]/F"   ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsLDisc")) tree->Branch((name+"Jet_DeepFlavourCvsLDisc").c_str(),  Jet_DeepFlavourCvsLDisc , (name+"Jet_DeepFlavourCvsLDisc["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsBDisc")) tree->Branch((name+"Jet_DeepFlavourCvsBDisc").c_str(),  Jet_DeepFlavourCvsBDisc , (name+"Jet_DeepFlavourCvsBDisc["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourB")) tree->Branch((name+"Jet_DeepFlavourB"       ).c_str(),  Jet_DeepFlavourB        , (name+"Jet_DeepFlavourB["+name+"nJet]/F"       ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourBB")) tree->Branch((name+"Jet_DeepFlavourBB"      ).c_str(),  Jet_DeepFlavourBB       , (name+"Jet_DeepFlavourBB["+name+"nJet]/F"      ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourLEPB")) tree->Branch((name+"Jet_DeepFlavourLEPB"    ).c_str(),  Jet_DeepFlavourLEPB     , (name+"Jet_DeepFlavourLEPB["+name+"nJet]/F"    ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourC")) tree->Branch((name+"Jet_DeepFlavourC"       ).c_str(),  Jet_DeepFlavourC        , (name+"Jet_DeepFlavourC["+name+"nJet]/F"       ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourUDS")) tree->Branch((name+"Jet_DeepFlavourUDS"     ).c_str(),  Jet_DeepFlavourUDS      , (name+"Jet_DeepFlavourUDS["+name+"nJet]/F"     ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourG")) tree->Branch((name+"Jet_DeepFlavourG"       ).c_str(),  Jet_DeepFlavourG        , (name+"Jet_DeepFlavourG["+name+"nJet]/F"       ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourBDiscN")) tree->Branch((name+"Jet_DeepFlavourBDiscN"   ).c_str(), Jet_DeepFlavourBDiscN   , (name+"Jet_DeepFlavourBDiscN["+name+"nJet]/F"   ).c_str());
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsLDiscN")) tree->Branch((name+"Jet_DeepFlavourCvsLDiscN").c_str(), Jet_DeepFlavourCvsLDiscN, (name+"Jet_DeepFlavourCvsLDiscN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsBDiscN")) tree->Branch((name+"Jet_DeepFlavourCvsBDiscN").c_str(), Jet_DeepFlavourCvsBDiscN, (name+"Jet_DeepFlavourCvsBDiscN["+name+"nJet]/F").c_str());


      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedBDisc")) tree->Branch((name+"Jet_DeepFlavourPrunedBDisc"   ).c_str(),  Jet_DeepFlavourPrunedBDisc    , (name+"Jet_DeepFlavourPrunedBDisc["+name+"nJet]/F"   ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsLDisc")) tree->Branch((name+"Jet_DeepFlavourPrunedCvsLDisc").c_str(),  Jet_DeepFlavourPrunedCvsLDisc , (name+"Jet_DeepFlavourPrunedCvsLDisc["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsBDisc")) tree->Branch((name+"Jet_DeepFlavourPrunedCvsBDisc").c_str(),  Jet_DeepFlavourPrunedCvsBDisc , (name+"Jet_DeepFlavourPrunedCvsBDisc["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedB")) tree->Branch((name+"Jet_DeepFlavourPrunedB"       ).c_str(),  Jet_DeepFlavourPrunedB        , (name+"Jet_DeepFlavourPrunedB["+name+"nJet]/F"       ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedBB")) tree->Branch((name+"Jet_DeepFlavourPrunedBB"      ).c_str(),  Jet_DeepFlavourPrunedBB       , (name+"Jet_DeepFlavourPrunedBB["+name+"nJet]/F"      ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedLEPB")) tree->Branch((name+"Jet_DeepFlavourPrunedLEPB"    ).c_str(),  Jet_DeepFlavourPrunedLEPB     , (name+"Jet_DeepFlavourPrunedLEPB["+name+"nJet]/F"    ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedC")) tree->Branch((name+"Jet_DeepFlavourPrunedC"       ).c_str(),  Jet_DeepFlavourPrunedC        , (name+"Jet_DeepFlavourPrunedC["+name+"nJet]/F"       ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedUDS")) tree->Branch((name+"Jet_DeepFlavourPrunedUDS"     ).c_str(),  Jet_DeepFlavourPrunedUDS      , (name+"Jet_DeepFlavourPrunedUDS["+name+"nJet]/F"     ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedG")) tree->Branch((name+"Jet_DeepFlavourPrunedG"       ).c_str(),  Jet_DeepFlavourPrunedG        , (name+"Jet_DeepFlavourPrunedG["+name+"nJet]/F"       ).c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedBDiscN")) tree->Branch((name+"Jet_DeepFlavourPrunedBDiscN"   ).c_str(), Jet_DeepFlavourPrunedBDiscN   , (name+"Jet_DeepFlavourPrunedBDiscN["+name+"nJet]/F"   ).c_str());
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsLDiscN")) tree->Branch((name+"Jet_DeepFlavourPrunedCvsLDiscN").c_str(), Jet_DeepFlavourPrunedCvsLDiscN, (name+"Jet_DeepFlavourPrunedCvsLDiscN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsBDiscN")) tree->Branch((name+"Jet_DeepFlavourPrunedCvsBDiscN").c_str(), Jet_DeepFlavourPrunedCvsBDiscN, (name+"Jet_DeepFlavourPrunedCvsBDiscN["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_DeepCSVBDisc")) tree->Branch((name+"Jet_DeepCSVBDisc"      ).c_str(), Jet_DeepCSVBDisc     , (name+"Jet_DeepCSVBDisc["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVBDiscN")) tree->Branch((name+"Jet_DeepCSVBDiscN"     ).c_str(), Jet_DeepCSVBDiscN    , (name+"Jet_DeepCSVBDiscN["+name+"nJet]/F").c_str()    );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVBDiscP")) tree->Branch((name+"Jet_DeepCSVBDiscP"     ).c_str(), Jet_DeepCSVBDiscP    , (name+"Jet_DeepCSVBDiscP["+name+"nJet]/F").c_str()    );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsLDisc")) tree->Branch((name+"Jet_DeepCSVCvsLDisc"   ).c_str(), Jet_DeepCSVCvsLDisc  , (name+"Jet_DeepCSVCvsLDisc["+name+"nJet]/F").c_str()  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsLDiscN")) tree->Branch((name+"Jet_DeepCSVCvsLDiscN"  ).c_str(), Jet_DeepCSVCvsLDiscN , (name+"Jet_DeepCSVCvsLDiscN["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsLDiscP")) tree->Branch((name+"Jet_DeepCSVCvsLDiscP"  ).c_str(), Jet_DeepCSVCvsLDiscP , (name+"Jet_DeepCSVCvsLDiscP["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsBDisc")) tree->Branch((name+"Jet_DeepCSVCvsBDisc"   ).c_str(), Jet_DeepCSVCvsBDisc  , (name+"Jet_DeepCSVCvsBDisc["+name+"nJet]/F").c_str()  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsBDiscN")) tree->Branch((name+"Jet_DeepCSVCvsBDiscN"  ).c_str(), Jet_DeepCSVCvsBDiscN , (name+"Jet_DeepCSVCvsBDiscN["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsBDiscP")) tree->Branch((name+"Jet_DeepCSVCvsBDiscP"  ).c_str(), Jet_DeepCSVCvsBDiscP , (name+"Jet_DeepCSVCvsBDiscP["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVb")) tree->Branch((name+"Jet_DeepCSVb"  ).c_str(),         Jet_DeepCSVb         , (name+"Jet_DeepCSVb["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVc")) tree->Branch((name+"Jet_DeepCSVc"  ).c_str(),         Jet_DeepCSVc         , (name+"Jet_DeepCSVc["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVl")) tree->Branch((name+"Jet_DeepCSVl"  ).c_str(),         Jet_DeepCSVl         , (name+"Jet_DeepCSVl["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbb")) tree->Branch((name+"Jet_DeepCSVbb" ).c_str(),         Jet_DeepCSVbb        , (name+"Jet_DeepCSVbb["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVcc")) tree->Branch((name+"Jet_DeepCSVcc" ).c_str(),         Jet_DeepCSVcc        , (name+"Jet_DeepCSVcc["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbN")) tree->Branch((name+"Jet_DeepCSVbN" ).c_str(),         Jet_DeepCSVbN        , (name+"Jet_DeepCSVbN["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVcN")) tree->Branch((name+"Jet_DeepCSVcN" ).c_str(),         Jet_DeepCSVcN        , (name+"Jet_DeepCSVcN["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVlN")) tree->Branch((name+"Jet_DeepCSVlN" ).c_str(),         Jet_DeepCSVlN        , (name+"Jet_DeepCSVlN["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbbN")) tree->Branch((name+"Jet_DeepCSVbbN").c_str(),         Jet_DeepCSVbbN       , (name+"Jet_DeepCSVbbN["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVccN")) tree->Branch((name+"Jet_DeepCSVccN").c_str(),         Jet_DeepCSVccN       , (name+"Jet_DeepCSVccN["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbP")) tree->Branch((name+"Jet_DeepCSVbP" ).c_str(),         Jet_DeepCSVbP        , (name+"Jet_DeepCSVbP["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVcP")) tree->Branch((name+"Jet_DeepCSVcP" ).c_str(),         Jet_DeepCSVcP        , (name+"Jet_DeepCSVcP["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVlP")) tree->Branch((name+"Jet_DeepCSVlP" ).c_str(),         Jet_DeepCSVlP        , (name+"Jet_DeepCSVlP["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbbP")) tree->Branch((name+"Jet_DeepCSVbbP").c_str(),         Jet_DeepCSVbbP       , (name+"Jet_DeepCSVbbP["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVccP")) tree->Branch((name+"Jet_DeepCSVccP").c_str(),         Jet_DeepCSVccP       , (name+"Jet_DeepCSVccP["+name+"nJet]/F").c_str()       );

      //--------------------------------------
      // pf electron information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstSE")) tree->Branch((name+"Jet_nFirstSE").c_str(),    Jet_nFirstSE    ,(name+"Jet_nFirstSE["+name+"nJet]/I").c_str() );
      if(variableParser.isToBeStored(name+"Jet_nLastSE")) tree->Branch((name+"Jet_nLastSE").c_str(),     Jet_nLastSE     ,(name+"Jet_nLastSE["+name+"nJet]/I").c_str()  );

      if(variableParser.isToBeStored(name+"nPFElectron")) tree->Branch((name+"nPFElectron").c_str()         ,&nPFElectron        ,(name+"nPFElectron/I").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_IdxJet")) tree->Branch((name+"PFElectron_IdxJet").c_str()   ,PFElectron_IdxJet   ,(name+"PFElectron_IdxJet["+name+"nPFElectron]/I").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_pt")) tree->Branch((name+"PFElectron_pt").c_str()       ,PFElectron_pt       ,(name+"PFElectron_pt["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_eta")) tree->Branch((name+"PFElectron_eta").c_str()      ,PFElectron_eta      ,(name+"PFElectron_eta["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_phi")) tree->Branch((name+"PFElectron_phi").c_str()      ,PFElectron_phi      ,(name+"PFElectron_phi["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_ptrel")) tree->Branch((name+"PFElectron_ptrel").c_str()    ,PFElectron_ptrel    ,(name+"PFElectron_ptrel["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_deltaR")) tree->Branch((name+"PFElectron_deltaR").c_str()   ,PFElectron_deltaR   ,(name+"PFElectron_deltaR["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_ratio")) tree->Branch((name+"PFElectron_ratio").c_str()    ,PFElectron_ratio    ,(name+"PFElectron_ratio["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_ratioRel")) tree->Branch((name+"PFElectron_ratioRel").c_str() ,PFElectron_ratioRel ,(name+"PFElectron_ratioRel["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_IP")) tree->Branch((name+"PFElectron_IP").c_str()       ,PFElectron_IP       ,(name+"PFElectron_IP["+name+"nPFElectron]/F").c_str());
      if(variableParser.isToBeStored(name+"PFElectron_IP2D")) tree->Branch((name+"PFElectron_IP2D").c_str()     ,PFElectron_IP2D     ,(name+"PFElectron_IP2D["+name+"nPFElectron]/F").c_str());

      //--------------------------------------
      // pf muon information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstSM")) tree->Branch((name+"Jet_nFirstSM").c_str()        ,Jet_nFirstSM        ,(name+"Jet_nFirstSM["+name+"nJet]/I").c_str() );
      if(variableParser.isToBeStored(name+"Jet_nLastSM")) tree->Branch((name+"Jet_nLastSM").c_str()         ,Jet_nLastSM         ,(name+"Jet_nLastSM["+name+"nJet]/I").c_str()  );

      if(variableParser.isToBeStored(name+"nPFMuon")) tree->Branch((name+"nPFMuon").c_str()            ,&nPFMuon            ,(name+"nPFMuon/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_IdxJet")) tree->Branch((name+"PFMuon_IdxJet").c_str()      ,PFMuon_IdxJet       ,(name+"PFMuon_IdxJet["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_nMuHit")) tree->Branch((name+"PFMuon_nMuHit").c_str()      ,PFMuon_nMuHit       ,(name+"PFMuon_nMuHit["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_nTkHit")) tree->Branch((name+"PFMuon_nTkHit").c_str()      ,PFMuon_nTkHit       ,(name+"PFMuon_nTkHit["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_nPixHit")) tree->Branch((name+"PFMuon_nPixHit").c_str()     ,PFMuon_nPixHit      ,(name+"PFMuon_nPixHit["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_nOutHit")) tree->Branch((name+"PFMuon_nOutHit").c_str()     ,PFMuon_nOutHit      ,(name+"PFMuon_nOutHit["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_nTkLwM")) tree->Branch((name+"PFMuon_nTkLwM").c_str()      ,PFMuon_nTkLwM       ,(name+"PFMuon_nTkLwM["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_nPixLwM")) tree->Branch((name+"PFMuon_nPixLwM").c_str()     ,PFMuon_nPixLwM      ,(name+"PFMuon_nPixLwM["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_nMatched")) tree->Branch((name+"PFMuon_nMatched").c_str()    ,PFMuon_nMatched     ,(name+"PFMuon_nMatched["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_chi2")) tree->Branch((name+"PFMuon_chi2").c_str()        ,PFMuon_chi2         ,(name+"PFMuon_chi2["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_chi2Tk")) tree->Branch((name+"PFMuon_chi2Tk").c_str()      ,PFMuon_chi2Tk       ,(name+"PFMuon_chi2Tk["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_isGlobal")) tree->Branch((name+"PFMuon_isGlobal").c_str()    ,PFMuon_isGlobal     ,(name+"PFMuon_isGlobal["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_hist")) tree->Branch((name+"PFMuon_hist").c_str()        ,PFMuon_hist         ,(name+"PFMuon_hist["+name+"nPFMuon]/I").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_pt")) tree->Branch((name+"PFMuon_pt").c_str()          ,PFMuon_pt           ,(name+"PFMuon_pt["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_eta")) tree->Branch((name+"PFMuon_eta").c_str()         ,PFMuon_eta          ,(name+"PFMuon_eta["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_phi")) tree->Branch((name+"PFMuon_phi").c_str()         ,PFMuon_phi          ,(name+"PFMuon_phi["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_ptrel")) tree->Branch((name+"PFMuon_ptrel").c_str()       ,PFMuon_ptrel        ,(name+"PFMuon_ptrel["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_deltaR")) tree->Branch((name+"PFMuon_deltaR").c_str()      ,PFMuon_deltaR       ,(name+"PFMuon_deltaR["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_ratio")) tree->Branch((name+"PFMuon_ratio").c_str()       ,PFMuon_ratio        ,(name+"PFMuon_ratio["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_ratioRel")) tree->Branch((name+"PFMuon_ratioRel").c_str()    ,PFMuon_ratioRel     ,(name+"PFMuon_ratioRel["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_IP")) tree->Branch((name+"PFMuon_IP").c_str()          ,PFMuon_IP           ,(name+"PFMuon_IP["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_IP2D")) tree->Branch((name+"PFMuon_IP2D").c_str()        ,PFMuon_IP2D         ,(name+"PFMuon_IP2D["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_IPsig")) tree->Branch((name+"PFMuon_IPsig").c_str()       ,PFMuon_IPsig        ,(name+"PFMuon_IPsig["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_IP2Dsig")) tree->Branch((name+"PFMuon_IP2Dsig").c_str()     ,PFMuon_IP2Dsig      ,(name+"PFMuon_IP2Dsig["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_dz")) tree->Branch((name+"PFMuon_dz").c_str()          ,PFMuon_dz           ,(name+"PFMuon_dz["+name+"nPFMuon]/F").c_str());
      if(variableParser.isToBeStored(name+"PFMuon_GoodQuality")) tree->Branch((name+"PFMuon_GoodQuality").c_str() ,PFMuon_GoodQuality  ,(name+"PFMuon_GoodQuality["+name+"nPFMuon]/I").c_str());

      //--------------------------------------
      // secondary vertex information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstSV")) tree->Branch((name+"Jet_nFirstSV").c_str(),       Jet_nFirstSV    ,(name+"Jet_nFirstSV["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nLastSV")) tree->Branch((name+"Jet_nLastSV").c_str(),        Jet_nLastSV     ,(name+"Jet_nLastSV["+name+"nJet]/I").c_str());

      if(variableParser.isToBeStored(name+"nSV")) tree->Branch((name+"nSV").c_str()                ,&nSV               ,(name+"nSV/I").c_str());
      if(variableParser.isToBeStored(name+"SV_x")) tree->Branch((name+"SV_x").c_str()               ,SV_x               ,(name+"SV_x["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_y")) tree->Branch((name+"SV_y").c_str()               ,SV_y               ,(name+"SV_y["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_z")) tree->Branch((name+"SV_z").c_str()               ,SV_z               ,(name+"SV_z["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_ex")) tree->Branch((name+"SV_ex").c_str()              ,SV_ex              ,(name+"SV_ex["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_ey")) tree->Branch((name+"SV_ey").c_str()              ,SV_ey              ,(name+"SV_ey["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_ez")) tree->Branch((name+"SV_ez").c_str()              ,SV_ez              ,(name+"SV_ez["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_chi2")) tree->Branch((name+"SV_chi2").c_str()            ,SV_chi2            ,(name+"SV_chi2["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_ndf")) tree->Branch((name+"SV_ndf").c_str()             ,SV_ndf             ,(name+"SV_ndf["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_flight")) tree->Branch((name+"SV_flight").c_str()          ,SV_flight          ,(name+"SV_flight["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_flightErr")) tree->Branch((name+"SV_flightErr").c_str()       ,SV_flightErr       ,(name+"SV_flightErr["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_deltaR_jet")) tree->Branch((name+"SV_deltaR_jet").c_str()      ,SV_deltaR_jet      ,(name+"SV_deltaR_jet["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_deltaR_sum_jet")) tree->Branch((name+"SV_deltaR_sum_jet").c_str()  ,SV_deltaR_sum_jet  ,(name+"SV_deltaR_sum_jet["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_deltaR_sum_dir")) tree->Branch((name+"SV_deltaR_sum_dir").c_str()  ,SV_deltaR_sum_dir  ,(name+"SV_deltaR_sum_dir["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_vtx_pt")) tree->Branch((name+"SV_vtx_pt").c_str()          ,SV_vtx_pt          ,(name+"SV_vtx_pt["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_flight2D")) tree->Branch((name+"SV_flight2D").c_str()        ,SV_flight2D        ,(name+"SV_flight2D["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_flight2DErr")) tree->Branch((name+"SV_flight2DErr").c_str()     ,SV_flight2DErr     ,(name+"SV_flight2DErr["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_totCharge")) tree->Branch((name+"SV_totCharge").c_str()       ,SV_totCharge       ,(name+"SV_totCharge ["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_vtxDistJetAxis")) tree->Branch((name+"SV_vtxDistJetAxis").c_str()  ,SV_vtxDistJetAxis  ,(name+"SV_vtxDistJetAxis ["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_nTrk")) tree->Branch((name+"SV_nTrk").c_str()            ,SV_nTrk            ,(name+"SV_nTrk["+name+"nSV]/I").c_str());
      if(variableParser.isToBeStored(name+"SV_mass")) tree->Branch((name+"SV_mass").c_str()            ,SV_mass            ,(name+"SV_mass["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_vtx_eta")) tree->Branch((name+"SV_vtx_eta").c_str()         ,SV_vtx_eta         ,(name+"SV_vtx_eta["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_vtx_phi")) tree->Branch((name+"SV_vtx_phi").c_str()         ,SV_vtx_phi         ,(name+"SV_vtx_phi["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_EnergyRatio")) tree->Branch((name+"SV_EnergyRatio").c_str()     ,SV_EnergyRatio     ,(name+"SV_EnergyRatio ["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_dir_x")) tree->Branch((name+"SV_dir_x").c_str()           ,SV_dir_x           ,(name+"SV_dir_x ["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_dir_y")) tree->Branch((name+"SV_dir_y").c_str()           ,SV_dir_y           ,(name+"SV_dir_y ["+name+"nSV]/F").c_str());
      if(variableParser.isToBeStored(name+"SV_dir_z")) tree->Branch((name+"SV_dir_z").c_str()           ,SV_dir_z           ,(name+"SV_dir_z ["+name+"nSV]/F").c_str());

      //--------------------------------------
      // track information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrack")) tree->Branch((name+"Jet_nFirstTrack").c_str(),  Jet_nFirstTrack ,(name+"Jet_nFirstTrack["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nLastTrack")) tree->Branch((name+"Jet_nLastTrack").c_str(),   Jet_nLastTrack  ,(name+"Jet_nLastTrack["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"nTrack")) tree->Branch((name+"nTrack").c_str()           ,&nTrack          ,(name+"nTrack/I").c_str());
      if(variableParser.isToBeStored(name+"Track_dxy")) tree->Branch((name+"Track_dxy").c_str()        ,Track_dxy             ,(name+"Track_dxy["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_dz")) tree->Branch((name+"Track_dz").c_str()         ,Track_dz         ,(name+"Track_dz["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_dxyError")) tree->Branch((name+"Track_dxyError").c_str()   ,Track_dxyError   ,(name+"Track_dxyError["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_dzError")) tree->Branch((name+"Track_dzError").c_str()    ,Track_dzError    ,(name+"Track_dzError["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_sign2D")) tree->Branch((name+"Track_sign2D").c_str()     ,Track_sign2D    ,(name+"Track_sign2D["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_sign3D")) tree->Branch((name+"Track_sign3D").c_str()     ,Track_sign3D    ,(name+"Track_sign3D["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_length")) tree->Branch((name+"Track_length").c_str()     ,Track_length     ,(name+"Track_length["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_dist")) tree->Branch((name+"Track_dist").c_str()       ,Track_dist            ,(name+"Track_dist["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_IP2D")) tree->Branch((name+"Track_IP2D").c_str()       ,Track_IP2D            ,(name+"Track_IP2D["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_IP2Dsig")) tree->Branch((name+"Track_IP2Dsig").c_str()    ,Track_IP2Dsig    ,(name+"Track_IP2Dsig["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_IP2Derr")) tree->Branch((name+"Track_IP2Derr").c_str()    ,Track_IP2Derr    ,(name+"Track_IP2Derr["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_IP")) tree->Branch((name+"Track_IP").c_str()         ,Track_IP         ,(name+"Track_IP["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_IPsig")) tree->Branch((name+"Track_IPsig").c_str()      ,Track_IPsig      ,(name+"Track_IPsig["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_IPerr")) tree->Branch((name+"Track_IPerr").c_str()      ,Track_IPerr      ,(name+"Track_IPerr["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_Proba")) tree->Branch((name+"Track_Proba").c_str()      ,Track_Proba      ,(name+"Track_Proba["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_p")) tree->Branch((name+"Track_p").c_str()          ,Track_p          ,(name+"Track_p["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_pt")) tree->Branch((name+"Track_pt").c_str()         ,Track_pt         ,(name+"Track_pt["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_eta")) tree->Branch((name+"Track_eta").c_str()        ,Track_eta             ,(name+"Track_eta["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_phi")) tree->Branch((name+"Track_phi").c_str()        ,Track_phi             ,(name+"Track_phi["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_chi2")) tree->Branch((name+"Track_chi2").c_str()       ,Track_chi2            ,(name+"Track_chi2["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_charge")) tree->Branch((name+"Track_charge").c_str()     ,Track_charge     ,(name+"Track_charge["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_history")) tree->Branch((name+"Track_history").c_str()    ,Track_history    ,(name+"Track_history["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitStrip")) tree->Branch((name+"Track_nHitStrip").c_str()  ,Track_nHitStrip  ,(name+"Track_nHitStrip["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitPixel")) tree->Branch((name+"Track_nHitPixel").c_str()  ,Track_nHitPixel  ,(name+"Track_nHitPixel["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitAll")) tree->Branch((name+"Track_nHitAll").c_str()    ,Track_nHitAll    ,(name+"Track_nHitAll["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitTIB")) tree->Branch((name+"Track_nHitTIB").c_str()    ,Track_nHitTIB    ,(name+"Track_nHitTIB["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitTID")) tree->Branch((name+"Track_nHitTID").c_str()    ,Track_nHitTID    ,(name+"Track_nHitTID["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitTOB")) tree->Branch((name+"Track_nHitTOB").c_str()    ,Track_nHitTOB    ,(name+"Track_nHitTOB["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitTEC")) tree->Branch((name+"Track_nHitTEC").c_str()    ,Track_nHitTEC    ,(name+"Track_nHitTEC["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitPXB")) tree->Branch((name+"Track_nHitPXB").c_str()    ,Track_nHitPXB    ,(name+"Track_nHitPXB["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_nHitPXF")) tree->Branch((name+"Track_nHitPXF").c_str()    ,Track_nHitPXF    ,(name+"Track_nHitPXF["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_isHitL1")) tree->Branch((name+"Track_isHitL1").c_str()    ,Track_isHitL1    ,(name+"Track_isHitL1["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_PV")) tree->Branch((name+"Track_PV").c_str()         ,Track_PV         ,(name+"Track_PV["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_SV")) tree->Branch((name+"Track_SV").c_str()         ,Track_SV         ,(name+"Track_SV["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_PVweight")) tree->Branch((name+"Track_PVweight").c_str()   ,Track_PVweight   ,(name+"Track_PVweight["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_SVweight")) tree->Branch((name+"Track_SVweight").c_str()   ,Track_SVweight   ,(name+"Track_SVweight["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_isfromSV")) tree->Branch((name+"Track_isfromSV").c_str()   ,Track_isfromSV   ,(name+"Track_isfromSV["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_isfromV0")) tree->Branch((name+"Track_isfromV0").c_str()   ,Track_isfromV0   ,(name+"Track_isfromV0["+name+"nTrack]/I").c_str());
      if(variableParser.isToBeStored(name+"Track_category")) tree->Branch((name+"Track_category").c_str()   ,Track_category   ,(name+"Track_category["+name+"nTrack]/I").c_str());

      //--------------------------------------
      // track truth information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Track_TPAssociationQuality")) tree->Branch((name+"Track_TPAssociationQuality").c_str()          ,Track_TPAssociationQuality             ,(name+"Track_TPAssociationQuality["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_idxMatchedTP")) tree->Branch((name+"Track_idxMatchedTP").c_str()                  ,Track_idxMatchedTP                     ,(name+"Track_idxMatchedTP["+name+"nTrack]/I").c_str());

      if(variableParser.isToBeStored(name+"Jet_nFirstTrackTruth")) tree->Branch((name+"Jet_nFirstTrackTruth").c_str()                ,Jet_nFirstTrackTruth                   ,(name+"Jet_nFirstTrackTruth["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nLastTrackTruth")) tree->Branch((name+"Jet_nLastTrackTruth").c_str()                 ,Jet_nLastTrackTruth                    ,(name+"Jet_nLastTrackTruth["+name+"nJet]/I").c_str());

      if(variableParser.isToBeStored(name+"nTrackTruth")) tree->Branch((name+"nTrackTruth").c_str()                ,&nTrackTruth                           ,(name+"nTrackTruth/I").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_idxMatchedTrack")) tree->Branch((name+"TrackTruth_idxMatchedTrack").c_str()          ,TrackTruth_idxMatchedTrack             ,(name+"TrackTruth_idxMatchedTrack["+name+"nTrackTruth]/I").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_p")) tree->Branch((name+"TrackTruth_p").c_str()                        ,TrackTruth_p                           ,(name+"TrackTruth_p["+name+"nTrackTruth]/F").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_pt")) tree->Branch((name+"TrackTruth_pt").c_str()                       ,TrackTruth_pt                          ,(name+"TrackTruth_pt["+name+"nTrackTruth]/F").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_eta")) tree->Branch((name+"TrackTruth_eta").c_str()                      ,TrackTruth_eta                         ,(name+"TrackTruth_eta["+name+"nTrackTruth]/F").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_phi")) tree->Branch((name+"TrackTruth_phi").c_str()                      ,TrackTruth_phi                         ,(name+"TrackTruth_phi["+name+"nTrackTruth]/F").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_charge")) tree->Branch((name+"TrackTruth_charge").c_str()                   ,TrackTruth_charge                      ,(name+"TrackTruth_charge["+name+"nTrackTruth]/I").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_pdgid")) tree->Branch((name+"TrackTruth_pdgid").c_str()                    ,TrackTruth_pdgid                       ,(name+"TrackTruth_pdgid["+name+"nTrackTruth]/I").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_dxy")) tree->Branch((name+"TrackTruth_dxy").c_str()                      ,TrackTruth_dxy                         ,(name+"TrackTruth_dxy["+name+"nTrackTruth]/F").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_dz")) tree->Branch((name+"TrackTruth_dz").c_str()                       ,TrackTruth_dz                          ,(name+"TrackTruth_dz["+name+"nTrackTruth]/F").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_nHitAll")) tree->Branch((name+"TrackTruth_nHitAll").c_str()                  ,TrackTruth_nHitAll                     ,(name+"TrackTruth_nHitAll["+name+"nTrackTruth]/I").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_nHitPixel")) tree->Branch((name+"TrackTruth_nHitPixel").c_str()                ,TrackTruth_nHitPixel                   ,(name+"TrackTruth_nHitPixel["+name+"nTrackTruth]/I").c_str());
      if(variableParser.isToBeStored(name+"TrackTruth_nHitStrip")) tree->Branch((name+"TrackTruth_nHitStrip").c_str()                ,TrackTruth_nHitStrip                   ,(name+"TrackTruth_nHitStrip["+name+"nTrackTruth]/I").c_str());

      //--------------------------------------
      // Inclusive Track information for PtRel template
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkInc")) tree->Branch((name+"Jet_nFirstTrkInc").c_str(), Jet_nFirstTrkInc ,(name+"Jet_nFirstTrkInc["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nLastTrkInc")) tree->Branch((name+"Jet_nLastTrkInc").c_str(),  Jet_nLastTrkInc  ,(name+"Jet_nLastTrkInc["+name+"nJet]/I").c_str());

      if(variableParser.isToBeStored(name+"nTrkInc")) tree->Branch((name+"nTrkInc").c_str()      ,&nTrkInc     ,(name+"nTrkInc/I").c_str());
      if(variableParser.isToBeStored(name+"TrkInc_pt")) tree->Branch((name+"TrkInc_pt").c_str()    ,TrkInc_pt    ,(name+"TrkInc_pt["+name+"nTrkInc]/F").c_str());
      if(variableParser.isToBeStored(name+"TrkInc_eta")) tree->Branch((name+"TrkInc_eta").c_str()   ,TrkInc_eta   ,(name+"TrkInc_eta["+name+"nTrkInc]/F").c_str());
      if(variableParser.isToBeStored(name+"TrkInc_phi")) tree->Branch((name+"TrkInc_phi").c_str()   ,TrkInc_phi   ,(name+"TrkInc_phi["+name+"nTrkInc]/F").c_str());
      if(variableParser.isToBeStored(name+"TrkInc_ptrel")) tree->Branch((name+"TrkInc_ptrel").c_str() ,TrkInc_ptrel ,(name+"TrkInc_ptrel["+name+"nTrkInc]/F").c_str());
      if(variableParser.isToBeStored(name+"TrkInc_IPsig")) tree->Branch((name+"TrkInc_IPsig").c_str() ,TrkInc_IPsig ,(name+"TrkInc_IPsig["+name+"nTrkInc]/F").c_str());
      if(variableParser.isToBeStored(name+"TrkInc_IP")) tree->Branch((name+"TrkInc_IP").c_str()    ,TrkInc_IP    ,(name+"TrkInc_IP["+name+"nTrkInc]/F").c_str());

      //--------------------------------------
      // TagInfo TaggingVariables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkTagVar")) tree->Branch((name+"Jet_nFirstTrkTagVar").c_str()  ,Jet_nFirstTrkTagVar  ,(name+"Jet_nFirstTrkTagVar["+name+"nJet]/I").c_str() );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkTagVar")) tree->Branch((name+"Jet_nLastTrkTagVar").c_str()   ,Jet_nLastTrkTagVar   ,(name+"Jet_nLastTrack["+name+"nJet]/I").c_str()      );
      if(variableParser.isToBeStored(name+"Jet_nFirstSVTagVar")) tree->Branch((name+"Jet_nFirstSVTagVar").c_str()   ,Jet_nFirstSVTagVar   ,(name+"Jet_nFirstSVTagVar["+name+"nJet]/I").c_str()  );
      if(variableParser.isToBeStored(name+"Jet_nLastSVTagVar")) tree->Branch((name+"Jet_nLastSVTagVar").c_str()    ,Jet_nLastSVTagVar    ,(name+"Jet_nLastSVTagVar["+name+"nJet]/I").c_str()   );

      if(variableParser.isToBeStored(name+"TagVar_jetNTracks")) tree->Branch((name+"TagVar_jetNTracks").c_str()                   ,TagVar_jetNTracks                   ,(name+"TagVar_jetNTracks["+name+"nJet]/I").c_str()                   );
      if(variableParser.isToBeStored(name+"TagVar_jetNSecondaryVertices")) tree->Branch((name+"TagVar_jetNSecondaryVertices").c_str()        ,TagVar_jetNSecondaryVertices        ,(name+"TagVar_jetNSecondaryVertices["+name+"nJet]/I").c_str()        );
      //run 2 id: https://twiki.cern.ch/CMS/JetID13TeVRun2017
      if(variableParser.isToBeStored(name+"TagVar_chargedHadronEnergyFraction")) tree->Branch((name+"TagVar_chargedHadronEnergyFraction").c_str()  ,TagVar_chargedHadronEnergyFraction  ,(name+"TagVar_chargedHadronEnergyFraction["+name+"nJet]/F").c_str()  );
      if(variableParser.isToBeStored(name+"TagVar_neutralHadronEnergyFraction")) tree->Branch((name+"TagVar_neutralHadronEnergyFraction").c_str()  ,TagVar_neutralHadronEnergyFraction  ,(name+"TagVar_neutralHadronEnergyFraction["+name+"nJet]/F").c_str()  );
      if(variableParser.isToBeStored(name+"TagVar_photonEnergyFraction")) tree->Branch((name+"TagVar_photonEnergyFraction").c_str()         ,TagVar_photonEnergyFraction         ,(name+"TagVar_photonEnergyFraction["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"TagVar_electronEnergyFraction")) tree->Branch((name+"TagVar_electronEnergyFraction").c_str()       ,TagVar_electronEnergyFraction       ,(name+"TagVar_electronEnergyFraction["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"TagVar_muonEnergyFraction")) tree->Branch((name+"TagVar_muonEnergyFraction").c_str()           ,TagVar_muonEnergyFraction           ,(name+"TagVar_muonEnergyFraction["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"TagVar_chargedHadronMultiplicity")) tree->Branch((name+"TagVar_chargedHadronMultiplicity").c_str()    ,TagVar_chargedHadronMultiplicity    ,(name+"TagVar_chargedHadronMultiplicity["+name+"nJet]/I").c_str()    );
      if(variableParser.isToBeStored(name+"TagVar_neutralHadronMultiplicity")) tree->Branch((name+"TagVar_neutralHadronMultiplicity").c_str()    ,TagVar_neutralHadronMultiplicity    ,(name+"TagVar_neutralHadronMultiplicity["+name+"nJet]/I").c_str()    );
      if(variableParser.isToBeStored(name+"TagVar_photonMultiplicity")) tree->Branch((name+"TagVar_photonMultiplicity").c_str()           ,TagVar_photonMultiplicity           ,(name+"TagVar_photonMultiplicity["+name+"nJet]/I").c_str()           );
      if(variableParser.isToBeStored(name+"TagVar_electronMultiplicity")) tree->Branch((name+"TagVar_electronMultiplicity").c_str()         ,TagVar_electronMultiplicity         ,(name+"TagVar_electronMultiplicity["+name+"nJet]/I").c_str()         );
      if(variableParser.isToBeStored(name+"TagVar_muonMultiplicity")) tree->Branch((name+"TagVar_muonMultiplicity").c_str()             ,TagVar_muonMultiplicity             ,(name+"TagVar_muonMultiplicity["+name+"nJet]/I").c_str()             );
      if(variableParser.isToBeStored(name+"TagVar_neutralEmEnergyFraction")) tree->Branch((name+"TagVar_neutralEmEnergyFraction").c_str()      ,TagVar_neutralEmEnergyFraction      ,(name+"TagVar_neutralEmEnergyFraction["+name+"nJet]/F").c_str()      );
      if(variableParser.isToBeStored(name+"TagVar_numberOfDaughters")) tree->Branch((name+"TagVar_numberOfDaughters").c_str()            ,TagVar_numberOfDaughters            ,(name+"TagVar_numberOfDaughters["+name+"nJet]/I").c_str()            );
      if(variableParser.isToBeStored(name+"TagVar_chargedMultiplicity")) tree->Branch((name+"TagVar_chargedMultiplicity").c_str()          ,TagVar_chargedMultiplicity          ,(name+"TagVar_chargedMultiplicity["+name+"nJet]/I").c_str()          );
      if(variableParser.isToBeStored(name+"TagVar_chargedEmEnergyFraction")) tree->Branch((name+"TagVar_chargedEmEnergyFraction").c_str()      ,TagVar_chargedEmEnergyFraction      ,(name+"TagVar_chargedEmEnergyFraction["+name+"nJet]/F").c_str()      );
      if(variableParser.isToBeStored(name+"TagVar_neutralMultiplicity")) tree->Branch((name+"TagVar_neutralMultiplicity").c_str()          ,TagVar_neutralMultiplicity          ,(name+"TagVar_neutralMultiplicity["+name+"nJet]/I").c_str()          );

      if(variableParser.isToBeStored(name+"nTrkTagVar")) tree->Branch((name+"nTrkTagVar").c_str()               ,&nTrkTagVar              ,(name+"nTrkTagVar/I").c_str()                                  );
      if(variableParser.isToBeStored(name+"TagVar_trackMomentum")) tree->Branch((name+"TagVar_trackMomentum").c_str()     ,TagVar_trackMomentum     ,(name+"TagVar_trackMomentum["+name+"nTrkTagVar]/F").c_str()    );
      if(variableParser.isToBeStored(name+"TagVar_trackEta")) tree->Branch((name+"TagVar_trackEta").c_str()          ,TagVar_trackEta          ,(name+"TagVar_trackEta["+name+"nTrkTagVar]/F").c_str()         );
      if(variableParser.isToBeStored(name+"TagVar_trackPhi")) tree->Branch((name+"TagVar_trackPhi").c_str()          ,TagVar_trackPhi          ,(name+"TagVar_trackPhi["+name+"nTrkTagVar]/F").c_str()         );
      if(variableParser.isToBeStored(name+"TagVar_trackPtRel")) tree->Branch((name+"TagVar_trackPtRel").c_str()        ,TagVar_trackPtRel        ,(name+"TagVar_trackPtRel["+name+"nTrkTagVar]/F").c_str()       );
      if(variableParser.isToBeStored(name+"TagVar_trackPPar")) tree->Branch((name+"TagVar_trackPPar").c_str()         ,TagVar_trackPPar         ,(name+"TagVar_trackPPar["+name+"nTrkTagVar]/F").c_str()        );
      if(variableParser.isToBeStored(name+"TagVar_trackEtaRel")) tree->Branch((name+"TagVar_trackEtaRel").c_str()       ,TagVar_trackEtaRel       ,(name+"TagVar_trackEtaRel["+name+"nTrkTagVar]/F").c_str()      );
      if(variableParser.isToBeStored(name+"TagVar_trackDeltaR")) tree->Branch((name+"TagVar_trackDeltaR").c_str()       ,TagVar_trackDeltaR       ,(name+"TagVar_trackDeltaR["+name+"nTrkTagVar]/F").c_str()      );
      if(variableParser.isToBeStored(name+"TagVar_trackPtRatio")) tree->Branch((name+"TagVar_trackPtRatio").c_str()      ,TagVar_trackPtRatio      ,(name+"TagVar_trackPtRatio["+name+"nTrkTagVar]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVar_trackPParRatio")) tree->Branch((name+"TagVar_trackPParRatio").c_str()    ,TagVar_trackPParRatio    ,(name+"TagVar_trackPParRatio["+name+"nTrkTagVar]/F").c_str()   );
      if(variableParser.isToBeStored(name+"TagVar_trackSip2dVal")) tree->Branch((name+"TagVar_trackSip2dVal").c_str()     ,TagVar_trackSip2dVal     ,(name+"TagVar_trackSip2dVal["+name+"nTrkTagVar]/F").c_str()    );
      if(variableParser.isToBeStored(name+"TagVar_trackSip2dSig")) tree->Branch((name+"TagVar_trackSip2dSig").c_str()     ,TagVar_trackSip2dSig     ,(name+"TagVar_trackSip2dSig["+name+"nTrkTagVar]/F").c_str()    );
      if(variableParser.isToBeStored(name+"TagVar_trackSip3dVal")) tree->Branch((name+"TagVar_trackSip3dVal").c_str()     ,TagVar_trackSip3dVal     ,(name+"TagVar_trackSip3dVal["+name+"nTrkTagVar]/F").c_str()    );
      if(variableParser.isToBeStored(name+"TagVar_trackSip3dSig")) tree->Branch((name+"TagVar_trackSip3dSig").c_str()     ,TagVar_trackSip3dSig     ,(name+"TagVar_trackSip3dSig["+name+"nTrkTagVar]/F").c_str()    );
      if(variableParser.isToBeStored(name+"TagVar_trackDecayLenVal")) tree->Branch((name+"TagVar_trackDecayLenVal").c_str()  ,TagVar_trackDecayLenVal  ,(name+"TagVar_trackDecayLenVal["+name+"nTrkTagVar]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVar_trackDecayLenSig")) tree->Branch((name+"TagVar_trackDecayLenSig").c_str()  ,TagVar_trackDecayLenSig  ,(name+"TagVar_trackDecayLenSig["+name+"nTrkTagVar]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVar_trackJetDistVal")) tree->Branch((name+"TagVar_trackJetDistVal").c_str()   ,TagVar_trackJetDistVal   ,(name+"TagVar_trackJetDistVal["+name+"nTrkTagVar]/F").c_str()  );
      if(variableParser.isToBeStored(name+"TagVar_trackJetDistSig")) tree->Branch((name+"TagVar_trackJetDistSig").c_str()   ,TagVar_trackJetDistSig   ,(name+"TagVar_trackJetDistSig["+name+"nTrkTagVar]/F").c_str()  );
      if(variableParser.isToBeStored(name+"TagVar_trackChi2")) tree->Branch((name+"TagVar_trackChi2").c_str()         ,TagVar_trackChi2         ,(name+"TagVar_trackChi2["+name+"nTrkTagVar]/F").c_str()        );
      if(variableParser.isToBeStored(name+"TagVar_trackNTotalHits")) tree->Branch((name+"TagVar_trackNTotalHits").c_str()   ,TagVar_trackNTotalHits   ,(name+"TagVar_trackNTotalHits["+name+"nTrkTagVar]/F").c_str()  );
      if(variableParser.isToBeStored(name+"TagVar_trackNPixelHits")) tree->Branch((name+"TagVar_trackNPixelHits").c_str()   ,TagVar_trackNPixelHits   ,(name+"TagVar_trackNPixelHits["+name+"nTrkTagVar]/F").c_str()  );

      if(variableParser.isToBeStored(name+"nSVTagVar")) tree->Branch((name+"nSVTagVar").c_str()                       ,&nSVTagVar                      ,(name+"nSVTagVar/I").c_str()                                     );
      if(variableParser.isToBeStored(name+"TagVar_vertexMass")) tree->Branch((name+"TagVar_vertexMass").c_str()               ,TagVar_vertexMass               ,(name+"TagVar_vertexMass["+name+"nSVTagVar]/F").c_str()          );
      if(variableParser.isToBeStored(name+"TagVar_vertexNTracks")) tree->Branch((name+"TagVar_vertexNTracks").c_str()            ,TagVar_vertexNTracks            ,(name+"TagVar_vertexNTracks["+name+"nSVTagVar]/F").c_str()       );
      if(variableParser.isToBeStored(name+"TagVar_vertexJetDeltaR")) tree->Branch((name+"TagVar_vertexJetDeltaR").c_str()          ,TagVar_vertexJetDeltaR          ,(name+"TagVar_vertexJetDeltaR["+name+"nSVTagVar]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance2dVal")) tree->Branch((name+"TagVar_flightDistance2dVal").c_str()      ,TagVar_flightDistance2dVal      ,(name+"TagVar_flightDistance2dVal["+name+"nSVTagVar]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance2dSig")) tree->Branch((name+"TagVar_flightDistance2dSig").c_str()      ,TagVar_flightDistance2dSig      ,(name+"TagVar_flightDistance2dSig["+name+"nSVTagVar]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance3dVal")) tree->Branch((name+"TagVar_flightDistance3dVal").c_str()      ,TagVar_flightDistance3dVal      ,(name+"TagVar_flightDistance3dVal["+name+"nSVTagVar]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance3dSig")) tree->Branch((name+"TagVar_flightDistance3dSig").c_str()      ,TagVar_flightDistance3dSig      ,(name+"TagVar_flightDistance3dSig["+name+"nSVTagVar]/F").c_str() );

      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_Sip3dVal")) tree->Branch((name+"DeepFlavourInput_charged_Sip3dVal").c_str(), DeepFlavourInput_charged_Sip3dVal, (name+"DeepFlavourInput_charged_Sip3dVal["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_Sip3dSig")) tree->Branch((name+"DeepFlavourInput_charged_Sip3dSig").c_str(), DeepFlavourInput_charged_Sip3dSig, (name+"DeepFlavourInput_charged_Sip3dSig["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_quality")) tree->Branch((name+"DeepFlavourInput_charged_quality" ).c_str(), DeepFlavourInput_charged_quality , (name+"DeepFlavourInput_charged_quality[" +name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_chi2")) tree->Branch((name+"DeepFlavourInput_charged_chi2"    ).c_str(), DeepFlavourInput_charged_chi2    , (name+"DeepFlavourInput_charged_chi2["    +name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"DeepFlavourInput_neutral_drminsv")) tree->Branch((name+"DeepFlavourInput_neutral_drminsv").c_str(), DeepFlavourInput_neutral_drminsv, (name+"DeepFlavourInput_neutral_drminsv["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"DeepFlavourInput_neutral_hadFrac")) tree->Branch((name+"DeepFlavourInput_neutral_hadFrac").c_str(), DeepFlavourInput_neutral_hadFrac, (name+"DeepFlavourInput_neutral_hadFrac["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"DeepFlavourInput_neutral_ptrel")) tree->Branch((name+"DeepFlavourInput_neutral_ptrel  ").c_str(), DeepFlavourInput_neutral_ptrel  , (name+"DeepFlavourInput_neutral_ptrel["  +name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"DeepFlavourInput_sv_d3d")) tree->Branch((name+"DeepFlavourInput_sv_d3d"     ).c_str(), DeepFlavourInput_sv_d3d     , (name+"DeepFlavourInput_sv_d3d["     +name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"DeepFlavourInput_sv_d3dsig")) tree->Branch((name+"DeepFlavourInput_sv_d3dsig"  ).c_str(), DeepFlavourInput_sv_d3dsig  , (name+"DeepFlavourInput_sv_d3dsig["  +name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"DeepFlavourInput_sv_normchi2")) tree->Branch((name+"DeepFlavourInput_sv_normchi2").c_str(), DeepFlavourInput_sv_normchi2, (name+"DeepFlavourInput_sv_normchi2["+name+"nJet]/F").c_str());

      //--------------------------------------
      // CSV TaggingVariables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"TagVarCSV_trackJetPt")) tree->Branch((name+"TagVarCSV_trackJetPt").c_str()               ,TagVarCSV_trackJetPt               ,(name+"TagVarCSV_trackJetPt["+name+"nJet]/F").c_str()              );
      if(variableParser.isToBeStored(name+"TagVarCSV_jetNTracks")) tree->Branch((name+"TagVarCSV_jetNTracks").c_str()               ,TagVarCSV_jetNTracks               ,(name+"TagVarCSV_jetNTracks["+name+"nJet]/F").c_str()              );
      if(variableParser.isToBeStored(name+"TagVarCSV_jetNTracksEtaRel")) tree->Branch((name+"TagVarCSV_jetNTracksEtaRel").c_str()         ,TagVarCSV_jetNTracksEtaRel         ,(name+"TagVarCSV_jetNTracksEtaRel["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSumJetEtRatio")) tree->Branch((name+"TagVarCSV_trackSumJetEtRatio").c_str()       ,TagVarCSV_trackSumJetEtRatio       ,(name+"TagVarCSV_trackSumJetEtRatio["+name+"nJet]/F").c_str()      );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSumJetDeltaR")) tree->Branch((name+"TagVarCSV_trackSumJetDeltaR").c_str()        ,TagVarCSV_trackSumJetDeltaR        ,(name+"TagVarCSV_trackSumJetDeltaR["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dValAboveCharm")) tree->Branch((name+"TagVarCSV_trackSip2dValAboveCharm").c_str()  ,TagVarCSV_trackSip2dValAboveCharm  ,(name+"TagVarCSV_trackSip2dValAboveCharm["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dSigAboveCharm")) tree->Branch((name+"TagVarCSV_trackSip2dSigAboveCharm").c_str()  ,TagVarCSV_trackSip2dSigAboveCharm  ,(name+"TagVarCSV_trackSip2dSigAboveCharm["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dValAboveCharm")) tree->Branch((name+"TagVarCSV_trackSip3dValAboveCharm").c_str()  ,TagVarCSV_trackSip3dValAboveCharm  ,(name+"TagVarCSV_trackSip3dValAboveCharm["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dSigAboveCharm")) tree->Branch((name+"TagVarCSV_trackSip3dSigAboveCharm").c_str()  ,TagVarCSV_trackSip3dSigAboveCharm  ,(name+"TagVarCSV_trackSip3dSigAboveCharm["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexCategory")) tree->Branch((name+"TagVarCSV_vertexCategory").c_str()           ,TagVarCSV_vertexCategory           ,(name+"TagVarCSV_vertexCategory["+name+"nJet]/F").c_str()          );
      if(variableParser.isToBeStored(name+"TagVarCSV_jetNSecondaryVertices")) tree->Branch((name+"TagVarCSV_jetNSecondaryVertices").c_str()    ,TagVarCSV_jetNSecondaryVertices    ,(name+"TagVarCSV_jetNSecondaryVertices["+name+"nJet]/F").c_str()   );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexMass")) tree->Branch((name+"TagVarCSV_vertexMass").c_str()               ,TagVarCSV_vertexMass               ,(name+"TagVarCSV_vertexMass["+name+"nJet]/F").c_str()              );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexNTracks")) tree->Branch((name+"TagVarCSV_vertexNTracks").c_str()            ,TagVarCSV_vertexNTracks            ,(name+"TagVarCSV_vertexNTracks["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexEnergyRatio")) tree->Branch((name+"TagVarCSV_vertexEnergyRatio").c_str()        ,TagVarCSV_vertexEnergyRatio        ,(name+"TagVarCSV_vertexEnergyRatio["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexJetDeltaR")) tree->Branch((name+"TagVarCSV_vertexJetDeltaR").c_str()          ,TagVarCSV_vertexJetDeltaR          ,(name+"TagVarCSV_vertexJetDeltaR["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance2dVal")) tree->Branch((name+"TagVarCSV_flightDistance2dVal").c_str()      ,TagVarCSV_flightDistance2dVal      ,(name+"TagVarCSV_flightDistance2dVal["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance2dSig")) tree->Branch((name+"TagVarCSV_flightDistance2dSig").c_str()      ,TagVarCSV_flightDistance2dSig      ,(name+"TagVarCSV_flightDistance2dSig["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance3dVal")) tree->Branch((name+"TagVarCSV_flightDistance3dVal").c_str()      ,TagVarCSV_flightDistance3dVal      ,(name+"TagVarCSV_flightDistance3dVal["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance3dSig")) tree->Branch((name+"TagVarCSV_flightDistance3dSig").c_str()      ,TagVarCSV_flightDistance3dSig      ,(name+"TagVarCSV_flightDistance3dSig["+name+"nJet]/F").c_str()     );

      //--------------------------------------
      // CSV Tagging Track Variables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkTagVarCSV")) tree->Branch((name+"Jet_nFirstTrkTagVarCSV").c_str()        ,Jet_nFirstTrkTagVarCSV        ,(name+"Jet_nFirstTrkTagVarCSV["+name+"nJet]/I").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkTagVarCSV")) tree->Branch((name+"Jet_nLastTrkTagVarCSV").c_str()         ,Jet_nLastTrkTagVarCSV         ,(name+"Jet_nLastTrkTagVarCSV["+name+"nJet]/I").c_str()        );
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkEtaRelTagVarCSV")) tree->Branch((name+"Jet_nFirstTrkEtaRelTagVarCSV").c_str()  ,Jet_nFirstTrkEtaRelTagVarCSV  ,(name+"Jet_nFirstTrkEtaRelTagVarCSV["+name+"nJet]/I").c_str() );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkEtaRelTagVarCSV")) tree->Branch((name+"Jet_nLastTrkEtaRelTagVarCSV").c_str()   ,Jet_nLastTrkEtaRelTagVarCSV   ,(name+"Jet_nLastTrkEtaRelTagVarCSV["+name+"nJet]/I").c_str()      );

      if(variableParser.isToBeStored(name+"nTrkTagVarCSV")) tree->Branch((name+"nTrkTagVarCSV").c_str()               ,&nTrkTagVarCSV              ,(name+"nTrkTagVarCSV/I").c_str()                                      );
      if(variableParser.isToBeStored(name+"nTrkEtaRelTagVarCSV")) tree->Branch((name+"nTrkEtaRelTagVarCSV").c_str()         ,&nTrkEtaRelTagVarCSV        ,(name+"nTrkEtaRelTagVarCSV/I").c_str()                                );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackMomentum")) tree->Branch((name+"TagVarCSV_trackMomentum").c_str()     ,TagVarCSV_trackMomentum     ,(name+"TagVarCSV_trackMomentum["+name+"nTrkTagVarCSV]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackEta")) tree->Branch((name+"TagVarCSV_trackEta").c_str()          ,TagVarCSV_trackEta          ,(name+"TagVarCSV_trackEta["+name+"nTrkTagVarCSV]/F").c_str()          );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPhi")) tree->Branch((name+"TagVarCSV_trackPhi").c_str()          ,TagVarCSV_trackPhi          ,(name+"TagVarCSV_trackPhi["+name+"nTrkTagVarCSV]/F").c_str()          );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPtRel")) tree->Branch((name+"TagVarCSV_trackPtRel").c_str()        ,TagVarCSV_trackPtRel        ,(name+"TagVarCSV_trackPtRel["+name+"nTrkTagVarCSV]/F").c_str()        );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPPar")) tree->Branch((name+"TagVarCSV_trackPPar").c_str()         ,TagVarCSV_trackPPar         ,(name+"TagVarCSV_trackPPar["+name+"nTrkTagVarCSV]/F").c_str()         );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackDeltaR")) tree->Branch((name+"TagVarCSV_trackDeltaR").c_str()       ,TagVarCSV_trackDeltaR       ,(name+"TagVarCSV_trackDeltaR["+name+"nTrkTagVarCSV]/F").c_str()       );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPtRatio")) tree->Branch((name+"TagVarCSV_trackPtRatio").c_str()      ,TagVarCSV_trackPtRatio      ,(name+"TagVarCSV_trackPtRatio["+name+"nTrkTagVarCSV]/F").c_str()      );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPParRatio")) tree->Branch((name+"TagVarCSV_trackPParRatio").c_str()    ,TagVarCSV_trackPParRatio    ,(name+"TagVarCSV_trackPParRatio["+name+"nTrkTagVarCSV]/F").c_str()    );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dVal")) tree->Branch((name+"TagVarCSV_trackSip2dVal").c_str()     ,TagVarCSV_trackSip2dVal     ,(name+"TagVarCSV_trackSip2dVal["+name+"nTrkTagVarCSV]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dSig")) tree->Branch((name+"TagVarCSV_trackSip2dSig").c_str()     ,TagVarCSV_trackSip2dSig     ,(name+"TagVarCSV_trackSip2dSig["+name+"nTrkTagVarCSV]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dVal")) tree->Branch((name+"TagVarCSV_trackSip3dVal").c_str()     ,TagVarCSV_trackSip3dVal     ,(name+"TagVarCSV_trackSip3dVal["+name+"nTrkTagVarCSV]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dSig")) tree->Branch((name+"TagVarCSV_trackSip3dSig").c_str()     ,TagVarCSV_trackSip3dSig     ,(name+"TagVarCSV_trackSip3dSig["+name+"nTrkTagVarCSV]/F").c_str()     );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackDecayLenVal")) tree->Branch((name+"TagVarCSV_trackDecayLenVal").c_str()  ,TagVarCSV_trackDecayLenVal  ,(name+"TagVarCSV_trackDecayLenVal["+name+"nTrkTagVarCSV]/F").c_str()  );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackDecayLenSig")) tree->Branch((name+"TagVarCSV_trackDecayLenSig").c_str()  ,TagVarCSV_trackDecayLenSig  ,(name+"TagVarCSV_trackDecayLenSig["+name+"nTrkTagVarCSV]/F").c_str()  );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackJetDistVal")) tree->Branch((name+"TagVarCSV_trackJetDistVal").c_str()   ,TagVarCSV_trackJetDistVal   ,(name+"TagVarCSV_trackJetDistVal["+name+"nTrkTagVarCSV]/F").c_str()   );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackJetDistSig")) tree->Branch((name+"TagVarCSV_trackJetDistSig").c_str()   ,TagVarCSV_trackJetDistSig   ,(name+"TagVarCSV_trackJetDistSig["+name+"nTrkTagVarCSV]/F").c_str()   );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackEtaRel")) tree->Branch((name+"TagVarCSV_trackEtaRel").c_str()       ,TagVarCSV_trackEtaRel       ,(name+"TagVarCSV_trackEtaRel["+name+"nTrkEtaRelTagVarCSV]/F").c_str() );

      //--------------------------------------
      // CTag TaggingVariables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsB")) tree->Branch((name+"CTag_Jet_CvsB").c_str(),     CTag_Jet_CvsB     ,(name+"CTag_Jet_CvsB["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsBN")) tree->Branch((name+"CTag_Jet_CvsBN").c_str(),    CTag_Jet_CvsBN    ,(name+"CTag_Jet_CvsBN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsBP")) tree->Branch((name+"CTag_Jet_CvsBP").c_str(),    CTag_Jet_CvsBP    ,(name+"CTag_Jet_CvsBP["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsL")) tree->Branch((name+"CTag_Jet_CvsL").c_str(),     CTag_Jet_CvsL     ,(name+"CTag_Jet_CvsL["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsLN")) tree->Branch((name+"CTag_Jet_CvsLN").c_str(),    CTag_Jet_CvsLN    ,(name+"CTag_Jet_CvsLN["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsLP")) tree->Branch((name+"CTag_Jet_CvsLP").c_str(),    CTag_Jet_CvsLP    ,(name+"CTag_Jet_CvsLP["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_nFirstTrkCTagVar")) tree->Branch((name+"Jet_nFirstTrkCTagVar").c_str()        ,Jet_nFirstTrkCTagVar        ,(name+"Jet_nFirstTrkCTagVar["+name+"nJet]/I").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkCTagVar")) tree->Branch((name+"Jet_nLastTrkCTagVar").c_str()         ,Jet_nLastTrkCTagVar         ,(name+"Jet_nLastTrkCTagVar["+name+"nJet]/I").c_str()            );
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkEtaRelCTagVar")) tree->Branch((name+"Jet_nFirstTrkEtaRelCTagVar").c_str()        ,Jet_nFirstTrkEtaRelCTagVar        ,(name+"Jet_nFirstTrkEtaRelCTagVar["+name+"nJet]/I").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkEtaRelCTagVar")) tree->Branch((name+"Jet_nLastTrkEtaRelCTagVar").c_str()        ,Jet_nLastTrkEtaRelCTagVar        ,(name+"Jet_nLastTrkEtaRelCTagVar["+name+"nJet]/I").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_nFirstLepCTagVar")) tree->Branch((name+"Jet_nFirstLepCTagVar").c_str()        ,Jet_nFirstLepCTagVar        ,(name+"Jet_nFirstLepCTagVar["+name+"nJet]/I").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_nLastLepCTagVar")) tree->Branch((name+"Jet_nLastLepCTagVar").c_str()         ,Jet_nLastLepCTagVar         ,(name+"Jet_nLastLepCTagVar["+name+"nJet]/I").c_str()            );

      if(variableParser.isToBeStored(name+"CTag_jetNTracks")) tree->Branch((name+"CTag_jetNTracks").c_str()               ,CTag_jetNTracks               ,(name+"CTag_jetNTracks["+name+"nJet]/F").c_str()              );
      if(variableParser.isToBeStored(name+"CTag_jetNTracksEtaRel")) tree->Branch((name+"CTag_jetNTracksEtaRel").c_str()         ,CTag_jetNTracksEtaRel         ,(name+"CTag_jetNTracksEtaRel["+name+"nJet]/F").c_str()        );
      if(variableParser.isToBeStored(name+"CTag_trackSumJetEtRatio")) tree->Branch((name+"CTag_trackSumJetEtRatio").c_str()       ,CTag_trackSumJetEtRatio       ,(name+"CTag_trackSumJetEtRatio["+name+"nJet]/F").c_str()      );
      if(variableParser.isToBeStored(name+"CTag_trackSumJetDeltaR")) tree->Branch((name+"CTag_trackSumJetDeltaR").c_str()        ,CTag_trackSumJetDeltaR        ,(name+"CTag_trackSumJetDeltaR["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"CTag_trackSip2dSigAboveCharm")) tree->Branch((name+"CTag_trackSip2dSigAboveCharm").c_str()  ,CTag_trackSip2dSigAboveCharm  ,(name+"CTag_trackSip2dSigAboveCharm["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"CTag_trackSip3dSigAboveCharm")) tree->Branch((name+"CTag_trackSip3dSigAboveCharm").c_str()  ,CTag_trackSip3dSigAboveCharm  ,(name+"CTag_trackSip3dSigAboveCharm["+name+"nJet]/F").c_str() );
      if(variableParser.isToBeStored(name+"CTag_vertexCategory")) tree->Branch((name+"CTag_vertexCategory").c_str()           ,CTag_vertexCategory           ,(name+"CTag_vertexCategory["+name+"nJet]/F").c_str()          );
      if(variableParser.isToBeStored(name+"CTag_jetNSecondaryVertices")) tree->Branch((name+"CTag_jetNSecondaryVertices").c_str()    ,CTag_jetNSecondaryVertices    ,(name+"CTag_jetNSecondaryVertices["+name+"nJet]/F").c_str()   );
      if(variableParser.isToBeStored(name+"CTag_vertexMass")) tree->Branch((name+"CTag_vertexMass").c_str()               ,CTag_vertexMass               ,(name+"CTag_vertexMass["+name+"nJet]/F").c_str()              );
      if(variableParser.isToBeStored(name+"CTag_vertexNTracks")) tree->Branch((name+"CTag_vertexNTracks").c_str()            ,CTag_vertexNTracks            ,(name+"CTag_vertexNTracks["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"CTag_vertexEnergyRatio")) tree->Branch((name+"CTag_vertexEnergyRatio").c_str()        ,CTag_vertexEnergyRatio        ,(name+"CTag_vertexEnergyRatio["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"CTag_vertexJetDeltaR")) tree->Branch((name+"CTag_vertexJetDeltaR").c_str()          ,CTag_vertexJetDeltaR          ,(name+"CTag_vertexJetDeltaR["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"CTag_flightDistance2dSig")) tree->Branch((name+"CTag_flightDistance2dSig").c_str()      ,CTag_flightDistance2dSig      ,(name+"CTag_flightDistance2dSig["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"CTag_flightDistance3dSig")) tree->Branch((name+"CTag_flightDistance3dSig").c_str()      ,CTag_flightDistance3dSig      ,(name+"CTag_flightDistance3dSig["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"CTag_massVertexEnergyFraction")) tree->Branch((name+"CTag_massVertexEnergyFraction").c_str(),     CTag_massVertexEnergyFraction     ,(name+"CTag_massVertexEnergyFraction["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_vertexBoostOverSqrtJetPt")) tree->Branch((name+"CTag_vertexBoostOverSqrtJetPt").c_str(),     CTag_vertexBoostOverSqrtJetPt     ,(name+"CTag_vertexBoostOverSqrtJetPt["+name+"nJet]/F").c_str());

      if(variableParser.isToBeStored(name+"nTrkCTagVar")) tree->Branch((name+"nTrkCTagVar").c_str()               ,&nTrkCTagVar              ,(name+"nTrkCTagVar/I").c_str()                                      );
      if(variableParser.isToBeStored(name+"nTrkEtaRelCTagVar")) tree->Branch((name+"nTrkEtaRelCTagVar").c_str()         ,&nTrkEtaRelCTagVar        ,(name+"nTrkEtaRelCTagVar/I").c_str()                                );
      if(variableParser.isToBeStored(name+"CTag_trackPtRel")) tree->Branch((name+"CTag_trackPtRel").c_str()        ,CTag_trackPtRel        ,(name+"CTag_trackPtRel["+name+"nTrkCTagVar]/F").c_str()        );
      if(variableParser.isToBeStored(name+"CTag_trackPPar")) tree->Branch((name+"CTag_trackPPar").c_str()         ,CTag_trackPPar         ,(name+"CTag_trackPPar["+name+"nTrkCTagVar]/F").c_str()         );
      if(variableParser.isToBeStored(name+"CTag_trackDeltaR")) tree->Branch((name+"CTag_trackDeltaR").c_str()       ,CTag_trackDeltaR       ,(name+"CTag_trackDeltaR["+name+"nTrkCTagVar]/F").c_str()       );
      if(variableParser.isToBeStored(name+"CTag_trackPtRatio")) tree->Branch((name+"CTag_trackPtRatio").c_str()      ,CTag_trackPtRatio      ,(name+"CTag_trackPtRatio["+name+"nTrkCTagVar]/F").c_str()      );
      if(variableParser.isToBeStored(name+"CTag_trackPParRatio")) tree->Branch((name+"CTag_trackPParRatio").c_str()    ,CTag_trackPParRatio    ,(name+"CTag_trackPParRatio["+name+"nTrkCTagVar]/F").c_str()    );
      if(variableParser.isToBeStored(name+"CTag_trackSip2dSig")) tree->Branch((name+"CTag_trackSip2dSig").c_str()     ,CTag_trackSip2dSig     ,(name+"CTag_trackSip2dSig["+name+"nTrkCTagVar]/F").c_str()     );
      if(variableParser.isToBeStored(name+"CTag_trackSip3dSig")) tree->Branch((name+"CTag_trackSip3dSig").c_str()     ,CTag_trackSip3dSig     ,(name+"CTag_trackSip3dSig["+name+"nTrkCTagVar]/F").c_str()     );
      if(variableParser.isToBeStored(name+"CTag_trackDecayLenVal")) tree->Branch((name+"CTag_trackDecayLenVal").c_str()  ,CTag_trackDecayLenVal  ,(name+"CTag_trackDecayLenVal["+name+"nTrkCTagVar]/F").c_str()  );
      if(variableParser.isToBeStored(name+"CTag_trackJetDistVal")) tree->Branch((name+"CTag_trackJetDistVal").c_str()   ,CTag_trackJetDistVal   ,(name+"CTag_trackJetDistVal["+name+"nTrkCTagVar]/F").c_str()   );
      if(variableParser.isToBeStored(name+"CTag_trackEtaRel")) tree->Branch((name+"CTag_trackEtaRel").c_str()       ,CTag_trackEtaRel       ,(name+"CTag_trackEtaRel["+name+"nTrkEtaRelCTagVar]/F").c_str() );

      if(variableParser.isToBeStored(name+"CTag_vertexLeptonCategory")) tree->Branch((name+"CTag_vertexLeptonCategory").c_str()    ,CTag_vertexLeptonCategory ,(name+"CTag_vertexLeptonCategory["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_jetNLeptons")) tree->Branch((name+"CTag_jetNLeptons").c_str()        ,CTag_jetNLeptons               ,(name+"CTag_jetNLeptons["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"nLeptons")) tree->Branch((name+"nLeptons").c_str()               ,&nLeptons              ,(name+"nLeptons/I").c_str());
      if(variableParser.isToBeStored(name+"CTag_leptonPtRel")) tree->Branch((name+"CTag_leptonPtRel").c_str()         ,CTag_leptonPtRel         ,(name+"CTag_leptonPtRel["+name+"nLeptons]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_leptonSip3d")) tree->Branch((name+"CTag_leptonSip3d").c_str()         ,CTag_leptonSip3d         ,(name+"CTag_leptonSip3d["+name+"nLeptons]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_leptonDeltaR")) tree->Branch((name+"CTag_leptonDeltaR").c_str()         ,CTag_leptonDeltaR         ,(name+"CTag_leptonDeltaR["+name+"nLeptons]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_leptonRatioRel")) tree->Branch((name+"CTag_leptonRatioRel").c_str()         ,CTag_leptonRatioRel         ,(name+"CTag_leptonRatioRel["+name+"nLeptons]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_leptonEtaRel")) tree->Branch((name+"CTag_leptonEtaRel").c_str()         ,CTag_leptonEtaRel         ,(name+"CTag_leptonEtaRel["+name+"nLeptons]/F").c_str());
      if(variableParser.isToBeStored(name+"CTag_leptonRatio")) tree->Branch((name+"CTag_leptonRatio").c_str()         ,CTag_leptonRatio         ,(name+"CTag_leptonRatio["+name+"nLeptons]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_FatJetIdx")) tree->Branch((name+"Jet_FatJetIdx").c_str(),   Jet_FatJetIdx   ,(name+"Jet_FatJetIdx["+name+"nJet]/I").c_str());

      if(variableParser.isToBeStored(name+"Jet_ptSoftDrop")) tree->Branch((name+"Jet_ptSoftDrop").c_str(),    Jet_ptSoftDrop    ,(name+"Jet_ptSoftDrop["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_etaSoftDrop")) tree->Branch((name+"Jet_etaSoftDrop").c_str(),   Jet_etaSoftDrop   ,(name+"Jet_etaSoftDrop["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_phiSoftDrop")) tree->Branch((name+"Jet_phiSoftDrop").c_str(),   Jet_phiSoftDrop   ,(name+"Jet_phiSoftDrop["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_massSoftDrop")) tree->Branch((name+"Jet_massSoftDrop").c_str(),  Jet_massSoftDrop  ,(name+"Jet_massSoftDrop["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_jecF0SoftDrop")) tree->Branch((name+"Jet_jecF0SoftDrop").c_str(), Jet_jecF0SoftDrop ,(name+"Jet_jecF0SoftDrop["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_ptPruned")) tree->Branch((name+"Jet_ptPruned").c_str(),    Jet_ptPruned    ,(name+"Jet_ptPruned["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_etaPruned")) tree->Branch((name+"Jet_etaPruned").c_str(),   Jet_etaPruned   ,(name+"Jet_etaPruned["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_phiPruned")) tree->Branch((name+"Jet_phiPruned").c_str(),   Jet_phiPruned   ,(name+"Jet_phiPruned["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_massPruned")) tree->Branch((name+"Jet_massPruned").c_str(),  Jet_massPruned  ,(name+"Jet_massPruned["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_jecF0Pruned")) tree->Branch((name+"Jet_jecF0Pruned").c_str(), Jet_jecF0Pruned ,(name+"Jet_jecF0Pruned["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tau1")) tree->Branch((name+"Jet_tau1").c_str(),        Jet_tau1        ,(name+"Jet_tau1["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tau2")) tree->Branch((name+"Jet_tau2").c_str(),        Jet_tau2        ,(name+"Jet_tau2["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tauAxis1_px")) tree->Branch((name+"Jet_tauAxis1_px").c_str(), Jet_tauAxis1_px ,(name+"Jet_tauAxis1_px["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tauAxis1_py")) tree->Branch((name+"Jet_tauAxis1_py").c_str(), Jet_tauAxis1_py ,(name+"Jet_tauAxis1_py["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tauAxis1_pz")) tree->Branch((name+"Jet_tauAxis1_pz").c_str(), Jet_tauAxis1_pz ,(name+"Jet_tauAxis1_pz["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tauAxis2_px")) tree->Branch((name+"Jet_tauAxis2_px").c_str(), Jet_tauAxis2_px ,(name+"Jet_tauAxis2_px["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tauAxis2_py")) tree->Branch((name+"Jet_tauAxis2_py").c_str(), Jet_tauAxis2_py ,(name+"Jet_tauAxis2_py["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_tauAxis2_pz")) tree->Branch((name+"Jet_tauAxis2_pz").c_str(), Jet_tauAxis2_pz ,(name+"Jet_tauAxis2_pz["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_z_ratio")) tree->Branch((name+"Jet_z_ratio").c_str(),          Jet_z_ratio          ,(name+"Jet_z_ratio["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_nTracks_fat")) tree->Branch((name+"Jet_nTracks_fat").c_str(),      Jet_nTracks_fat      ,(name+"Jet_nTracks_fat["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"Jet_nSV_fat")) tree->Branch((name+"Jet_nSV_fat").c_str(),          Jet_nSV_fat          ,(name+"Jet_nSV_fat["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackEtaRel_0")) tree->Branch((name+"Jet_tau1_trackEtaRel_0").c_str(),           Jet_tau1_trackEtaRel_0           ,(name+"Jet_tau1_trackEtaRel_0["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackEtaRel_1")) tree->Branch((name+"Jet_tau1_trackEtaRel_1").c_str(),           Jet_tau1_trackEtaRel_1           ,(name+"Jet_tau1_trackEtaRel_1["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackEtaRel_2")) tree->Branch((name+"Jet_tau1_trackEtaRel_2").c_str(),           Jet_tau1_trackEtaRel_2           ,(name+"Jet_tau1_trackEtaRel_2["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackEtaRel_0")) tree->Branch((name+"Jet_tau2_trackEtaRel_0").c_str(),           Jet_tau2_trackEtaRel_0           ,(name+"Jet_tau2_trackEtaRel_0["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackEtaRel_1")) tree->Branch((name+"Jet_tau2_trackEtaRel_1").c_str(),           Jet_tau2_trackEtaRel_1           ,(name+"Jet_tau2_trackEtaRel_1["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackEtaRel_2")) tree->Branch((name+"Jet_tau2_trackEtaRel_2").c_str(),           Jet_tau2_trackEtaRel_2           ,(name+"Jet_tau2_trackEtaRel_2["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_tau1_nSecondaryVertices")) tree->Branch((name+"Jet_tau1_nSecondaryVertices").c_str(),      Jet_tau1_nSecondaryVertices      ,(name+"Jet_tau1_nSecondaryVertices["+name+"nJet]/F").c_str()      );
      if(variableParser.isToBeStored(name+"Jet_tau2_nSecondaryVertices")) tree->Branch((name+"Jet_tau2_nSecondaryVertices").c_str(),      Jet_tau2_nSecondaryVertices      ,(name+"Jet_tau2_nSecondaryVertices["+name+"nJet]/F").c_str()      );
      if(variableParser.isToBeStored(name+"Jet_tau1_flightDistance2dSig")) tree->Branch((name+"Jet_tau1_flightDistance2dSig").c_str(),     Jet_tau1_flightDistance2dSig     ,(name+"Jet_tau1_flightDistance2dSig["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"Jet_tau2_flightDistance2dSig")) tree->Branch((name+"Jet_tau2_flightDistance2dSig").c_str(),     Jet_tau2_flightDistance2dSig     ,(name+"Jet_tau2_flightDistance2dSig["+name+"nJet]/F").c_str()     );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexDeltaR")) tree->Branch((name+"Jet_tau1_vertexDeltaR").c_str(),            Jet_tau1_vertexDeltaR            ,(name+"Jet_tau1_vertexDeltaR["+name+"nJet]/F").c_str()            );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexDeltaR")) tree->Branch((name+"Jet_tau2_vertexDeltaR").c_str(),            Jet_tau2_vertexDeltaR            ,(name+"Jet_tau2_vertexDeltaR["+name+"nJet]/F").c_str()            );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexEnergyRatio")) tree->Branch((name+"Jet_tau1_vertexEnergyRatio").c_str(),       Jet_tau1_vertexEnergyRatio       ,(name+"Jet_tau1_vertexEnergyRatio["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexEnergyRatio")) tree->Branch((name+"Jet_tau2_vertexEnergyRatio").c_str(),       Jet_tau2_vertexEnergyRatio       ,(name+"Jet_tau2_vertexEnergyRatio["+name+"nJet]/F").c_str()       );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexMass")) tree->Branch((name+"Jet_tau1_vertexMass").c_str(),              Jet_tau1_vertexMass              ,(name+"Jet_tau1_vertexMass["+name+"nJet]/F").c_str()              );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexMass")) tree->Branch((name+"Jet_tau2_vertexMass").c_str(),              Jet_tau2_vertexMass              ,(name+"Jet_tau2_vertexMass["+name+"nJet]/F").c_str()              );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexMass_corrected")) tree->Branch((name+"Jet_tau1_vertexMass_corrected").c_str(),    Jet_tau1_vertexMass_corrected    ,(name+"Jet_tau1_vertexMass_corrected["+name+"nJet]/F").c_str()    );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexMass_corrected")) tree->Branch((name+"Jet_tau2_vertexMass_corrected").c_str(),    Jet_tau2_vertexMass_corrected    ,(name+"Jet_tau2_vertexMass_corrected["+name+"nJet]/F").c_str()    );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexNTracks")) tree->Branch((name+"Jet_tau1_vertexNTracks").c_str(),           Jet_tau1_vertexNTracks           ,(name+"Jet_tau1_vertexNTracks["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexNTracks")) tree->Branch((name+"Jet_tau2_vertexNTracks").c_str(),           Jet_tau2_vertexNTracks           ,(name+"Jet_tau2_vertexNTracks["+name+"nJet]/F").c_str()           );
      if(variableParser.isToBeStored(name+"Jet_DoubleSV")) tree->Branch((name+"Jet_DoubleSV").c_str(),         Jet_DoubleSV         ,(name+"Jet_DoubleSV["+name+"nJet]/F").c_str());
      if(variableParser.isToBeStored(name+"Jet_DeepDoubleB")) tree->Branch((name+"Jet_DeepDoubleB").c_str(),      Jet_DeepDoubleB      ,(name+"Jet_DeepDoubleB["+name+"nJet]/F").c_str());


      if(variableParser.isToBeStored(name+"Track_lengthTau")) tree->Branch((name+"Track_lengthTau").c_str()     ,Track_lengthTau       ,(name+"Track_lengthTau["+name+"nTrack]/F").c_str());
      if(variableParser.isToBeStored(name+"Track_distTau")) tree->Branch((name+"Track_distTau").c_str()       ,Track_distTau         ,(name+"Track_distTau["+name+"nTrack]/F").c_str());

      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_3")) tree->Branch((name+"Jet_trackSip3dSig_3").c_str(),                Jet_trackSip3dSig_3                ,(name+"Jet_trackSip3dSig_3["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_2")) tree->Branch((name+"Jet_trackSip3dSig_2").c_str(),                Jet_trackSip3dSig_2                ,(name+"Jet_trackSip3dSig_2["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_1")) tree->Branch((name+"Jet_trackSip3dSig_1").c_str(),                Jet_trackSip3dSig_1                ,(name+"Jet_trackSip3dSig_1["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_0")) tree->Branch((name+"Jet_trackSip3dSig_0").c_str(),                Jet_trackSip3dSig_0                ,(name+"Jet_trackSip3dSig_0["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveCharm_0")) tree->Branch((name+"Jet_trackSip2dSigAboveCharm_0").c_str(),      Jet_trackSip2dSigAboveCharm_0      ,(name+"Jet_trackSip2dSigAboveCharm_0["+name+"nJet]/F").c_str()             );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveCharm_1")) tree->Branch((name+"Jet_trackSip2dSigAboveCharm_1").c_str(),      Jet_trackSip2dSigAboveCharm_1      ,(name+"Jet_trackSip2dSigAboveCharm_1["+name+"nJet]/F").c_str()             );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveBottom_0")) tree->Branch((name+"Jet_trackSip2dSigAboveBottom_0").c_str(),     Jet_trackSip2dSigAboveBottom_0     ,(name+"Jet_trackSip2dSigAboveBottom_0["+name+"nJet]/F").c_str()            );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveBottom_1")) tree->Branch((name+"Jet_trackSip2dSigAboveBottom_1").c_str(),     Jet_trackSip2dSigAboveBottom_1     ,(name+"Jet_trackSip2dSigAboveBottom_1["+name+"nJet]/F").c_str()            );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackSip3dSig_0")) tree->Branch((name+"Jet_tau1_trackSip3dSig_0").c_str(),           Jet_tau1_trackSip3dSig_0           ,(name+"Jet_tau1_trackSip3dSig_0["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackSip3dSig_1")) tree->Branch((name+"Jet_tau1_trackSip3dSig_1").c_str(),           Jet_tau1_trackSip3dSig_1           ,(name+"Jet_tau1_trackSip3dSig_1["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackSip3dSig_0")) tree->Branch((name+"Jet_tau2_trackSip3dSig_0").c_str(),           Jet_tau2_trackSip3dSig_0           ,(name+"Jet_tau2_trackSip3dSig_0["+name+"nJet]/F").c_str()         );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackSip3dSig_1")) tree->Branch((name+"Jet_tau2_trackSip3dSig_1").c_str(),           Jet_tau2_trackSip3dSig_1           ,(name+"Jet_tau2_trackSip3dSig_1["+name+"nJet]/F").c_str()         );
    }

    void ReadBranches(TTree *tree, VariableParser variableParser, std::string name="") {
      if(name!="") name += ".";
      if(variableParser.isToBeStored(name+"nJet")) tree->SetBranchAddress((name+"nJet").c_str(),            &nJet           );
      if(variableParser.isToBeStored(name+"Jet_pt")) tree->SetBranchAddress((name+"Jet_pt").c_str(),          Jet_pt          );
      if(variableParser.isToBeStored(name+"Jet_uncorrpt")) tree->SetBranchAddress((name+"Jet_uncorrpt").c_str(),    Jet_uncorrpt    );
      if(variableParser.isToBeStored(name+"Jet_genpt")) tree->SetBranchAddress((name+"Jet_genpt").c_str(),       Jet_genpt       );
      if(variableParser.isToBeStored(name+"Jet_residual")) tree->SetBranchAddress((name+"Jet_residual").c_str(),    Jet_residual    );
      if(variableParser.isToBeStored(name+"Jet_area")) tree->SetBranchAddress((name+"Jet_area").c_str(),        Jet_area        );
      if(variableParser.isToBeStored(name+"Jet_jes")) tree->SetBranchAddress((name+"Jet_jes").c_str(),         Jet_jes         );
      if(variableParser.isToBeStored(name+"Jet_eta")) tree->SetBranchAddress((name+"Jet_eta").c_str(),         Jet_eta         );
      if(variableParser.isToBeStored(name+"Jet_phi")) tree->SetBranchAddress((name+"Jet_phi").c_str(),         Jet_phi         );
      if(variableParser.isToBeStored(name+"Jet_mass")) tree->SetBranchAddress((name+"Jet_mass").c_str(),        Jet_mass        );
      if(variableParser.isToBeStored(name+"Jet_ntracks")) tree->SetBranchAddress((name+"Jet_ntracks").c_str(),     Jet_ntracks     );
      if(variableParser.isToBeStored(name+"Jet_nseltracks")) tree->SetBranchAddress((name+"Jet_nseltracks").c_str(),  Jet_nseltracks  );
      if(variableParser.isToBeStored(name+"Jet_flavour")) tree->SetBranchAddress((name+"Jet_flavour").c_str(),     Jet_flavour     );
      if(variableParser.isToBeStored(name+"Jet_flavourCleaned")) tree->SetBranchAddress((name+"Jet_flavourCleaned").c_str(), Jet_flavourCleaned);
      if(variableParser.isToBeStored(name+"Jet_partonFlavour")) tree->SetBranchAddress((name+"Jet_partonFlavour").c_str(), Jet_partonFlavour);
      if(variableParser.isToBeStored(name+"Jet_hadronFlavour")) tree->SetBranchAddress((name+"Jet_hadronFlavour").c_str(), Jet_hadronFlavour);
      if(variableParser.isToBeStored(name+"Jet_nbHadrons")) tree->SetBranchAddress((name+"Jet_nbHadrons").c_str(),   Jet_nbHadrons   );
      if(variableParser.isToBeStored(name+"Jet_ncHadrons")) tree->SetBranchAddress((name+"Jet_ncHadrons").c_str(),   Jet_ncHadrons   );
      if(variableParser.isToBeStored(name+"Jet_Ip2N")) tree->SetBranchAddress((name+"Jet_Ip2N").c_str(),        Jet_Ip2N        );
      if(variableParser.isToBeStored(name+"Jet_Ip2P")) tree->SetBranchAddress((name+"Jet_Ip2P").c_str(),        Jet_Ip2P        );
      if(variableParser.isToBeStored(name+"Jet_Ip3N")) tree->SetBranchAddress((name+"Jet_Ip3N").c_str(),        Jet_Ip3N        );
      if(variableParser.isToBeStored(name+"Jet_Ip3P")) tree->SetBranchAddress((name+"Jet_Ip3P").c_str(),        Jet_Ip3P        );
      if(variableParser.isToBeStored(name+"Jet_ProbaN")) tree->SetBranchAddress((name+"Jet_ProbaN").c_str(),      Jet_ProbaN     );
      if(variableParser.isToBeStored(name+"Jet_ProbaP")) tree->SetBranchAddress((name+"Jet_ProbaP").c_str(),      Jet_ProbaP     );
      if(variableParser.isToBeStored(name+"Jet_Proba")) tree->SetBranchAddress((name+"Jet_Proba").c_str(),       Jet_Proba      );
      if(variableParser.isToBeStored(name+"Jet_BprobN")) tree->SetBranchAddress((name+"Jet_BprobN").c_str(),      Jet_BprobN     );
      if(variableParser.isToBeStored(name+"Jet_BprobP")) tree->SetBranchAddress((name+"Jet_BprobP").c_str(),      Jet_BprobP     );
      if(variableParser.isToBeStored(name+"Jet_Bprob")) tree->SetBranchAddress((name+"Jet_Bprob").c_str(),       Jet_Bprob      );
      if(variableParser.isToBeStored(name+"Jet_SvxN")) tree->SetBranchAddress((name+"Jet_SvxN").c_str(),        Jet_SvxN       );
      if(variableParser.isToBeStored(name+"Jet_Svx")) tree->SetBranchAddress((name+"Jet_Svx").c_str(),         Jet_Svx        );
      if(variableParser.isToBeStored(name+"Jet_SvxNHP")) tree->SetBranchAddress((name+"Jet_SvxNHP").c_str(),      Jet_SvxNHP     );
      if(variableParser.isToBeStored(name+"Jet_SvxHP")) tree->SetBranchAddress((name+"Jet_SvxHP").c_str(),       Jet_SvxHP      );
      if(variableParser.isToBeStored(name+"Jet_CombSvxN")) tree->SetBranchAddress((name+"Jet_CombSvxN").c_str(),    Jet_CombSvxN   );
      if(variableParser.isToBeStored(name+"Jet_CombSvxP")) tree->SetBranchAddress((name+"Jet_CombSvxP").c_str(),    Jet_CombSvxP   );
      if(variableParser.isToBeStored(name+"Jet_CombSvx")) tree->SetBranchAddress((name+"Jet_CombSvx").c_str(),     Jet_CombSvx    );

      if(variableParser.isToBeStored(name+"Jet_CombIVF")) tree->SetBranchAddress((name+"Jet_CombIVF").c_str(),     Jet_CombIVF     );
      if(variableParser.isToBeStored(name+"Jet_CombIVF_P")) tree->SetBranchAddress((name+"Jet_CombIVF_P").c_str(),   Jet_CombIVF_P   );
      if(variableParser.isToBeStored(name+"Jet_CombIVF_N")) tree->SetBranchAddress((name+"Jet_CombIVF_N").c_str(),   Jet_CombIVF_N   );

      if(variableParser.isToBeStored(name+"Jet_SoftMuN")) tree->SetBranchAddress((name+"Jet_SoftMuN").c_str(),     Jet_SoftMuN     );
      if(variableParser.isToBeStored(name+"Jet_SoftMuP")) tree->SetBranchAddress((name+"Jet_SoftMuP").c_str(),     Jet_SoftMuP     );
      if(variableParser.isToBeStored(name+"Jet_SoftMu")) tree->SetBranchAddress((name+"Jet_SoftMu").c_str(),      Jet_SoftMu      );

      if(variableParser.isToBeStored(name+"Jet_SoftElN")) tree->SetBranchAddress((name+"Jet_SoftElN").c_str(),     Jet_SoftElN     );
      if(variableParser.isToBeStored(name+"Jet_SoftElP")) tree->SetBranchAddress((name+"Jet_SoftElP").c_str(),     Jet_SoftElP     );
      if(variableParser.isToBeStored(name+"Jet_SoftEl")) tree->SetBranchAddress((name+"Jet_SoftEl").c_str(),      Jet_SoftEl      );

      if(variableParser.isToBeStored(name+"Jet_cMVA")) tree->SetBranchAddress((name+"Jet_cMVA").c_str(),    Jet_cMVA    );
      if(variableParser.isToBeStored(name+"Jet_cMVAv2")) tree->SetBranchAddress((name+"Jet_cMVAv2").c_str(),    Jet_cMVAv2    );
      if(variableParser.isToBeStored(name+"Jet_cMVAv2N")) tree->SetBranchAddress((name+"Jet_cMVAv2N").c_str(),    Jet_cMVAv2N    );
      if(variableParser.isToBeStored(name+"Jet_cMVAv2P")) tree->SetBranchAddress((name+"Jet_cMVAv2P").c_str(),    Jet_cMVAv2P    );

      if(variableParser.isToBeStored(name+"Jet_hist1")) tree->SetBranchAddress((name+"Jet_hist1").c_str(),       Jet_hist1       );
      if(variableParser.isToBeStored(name+"Jet_hist2")) tree->SetBranchAddress((name+"Jet_hist2").c_str(),       Jet_hist2       );
      if(variableParser.isToBeStored(name+"Jet_hist3")) tree->SetBranchAddress((name+"Jet_hist3").c_str(),       Jet_hist3       );
      if(variableParser.isToBeStored(name+"Jet_histJet")) tree->SetBranchAddress((name+"Jet_histJet").c_str(),     Jet_histJet     );
      if(variableParser.isToBeStored(name+"Jet_histSvx")) tree->SetBranchAddress((name+"Jet_histSvx").c_str(),     Jet_histSvx     );

      if(variableParser.isToBeStored(name+"Jet_SV_multi")) tree->SetBranchAddress((name+"Jet_SV_multi").c_str(),    Jet_SV_multi      );
      if(variableParser.isToBeStored(name+"Jet_nSM")) tree->SetBranchAddress((name+"Jet_nSM").c_str(),         Jet_nSM               );
      if(variableParser.isToBeStored(name+"Jet_nSE")) tree->SetBranchAddress((name+"Jet_nSE").c_str(),         Jet_nSE               );

      if(variableParser.isToBeStored(name+"Jet_looseID")) tree->SetBranchAddress((name+"Jet_looseID").c_str(),      Jet_looseID  );
      if(variableParser.isToBeStored(name+"Jet_tightID")) tree->SetBranchAddress((name+"Jet_tightID").c_str(),      Jet_tightID  );

      //test variables for AK4 jets: clean up in the future
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSig_AboveBottom_0")) tree->SetBranchAddress((name+"Jet_trackSip2dSig_AboveBottom_0").c_str(),     Jet_trackSip2dSig_AboveBottom_0                 );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSig_AboveBottom_1")) tree->SetBranchAddress((name+"Jet_trackSip2dSig_AboveBottom_1").c_str(),     Jet_trackSip2dSig_AboveBottom_1                 );

      if(variableParser.isToBeStored(name+"Jet_DeepFlavourBDisc")) tree->SetBranchAddress((name+"Jet_DeepFlavourBDisc"   ).c_str(),  Jet_DeepFlavourBDisc     );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsLDisc")) tree->SetBranchAddress((name+"Jet_DeepFlavourCvsLDisc").c_str(),  Jet_DeepFlavourCvsLDisc  );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsBDisc")) tree->SetBranchAddress((name+"Jet_DeepFlavourCvsBDisc").c_str(),  Jet_DeepFlavourCvsBDisc  );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourB")) tree->SetBranchAddress((name+"Jet_DeepFlavourB"       ).c_str(),  Jet_DeepFlavourB         );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourBB")) tree->SetBranchAddress((name+"Jet_DeepFlavourBB"      ).c_str(),  Jet_DeepFlavourBB        );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourLEPB")) tree->SetBranchAddress((name+"Jet_DeepFlavourLEPB"    ).c_str(),  Jet_DeepFlavourLEPB      );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourC")) tree->SetBranchAddress((name+"Jet_DeepFlavourC"       ).c_str(),  Jet_DeepFlavourC         );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourUDS")) tree->SetBranchAddress((name+"Jet_DeepFlavourUDS"     ).c_str(),  Jet_DeepFlavourUDS       );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourG")) tree->SetBranchAddress((name+"Jet_DeepFlavourG"       ).c_str(),  Jet_DeepFlavourG         );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourBDiscN")) tree->SetBranchAddress((name+"Jet_DeepFlavourBDiscN"   ).c_str(), Jet_DeepFlavourBDiscN   );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsLDiscN")) tree->SetBranchAddress((name+"Jet_DeepFlavourCvsLDiscN").c_str(), Jet_DeepFlavourCvsLDiscN);
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourCvsBDiscN")) tree->SetBranchAddress((name+"Jet_DeepFlavourCvsBDiscN").c_str(), Jet_DeepFlavourCvsBDiscN);

      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedBDisc")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedBDisc"   ).c_str(),  Jet_DeepFlavourPrunedBDisc     );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsLDisc")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedCvsLDisc").c_str(),  Jet_DeepFlavourPrunedCvsLDisc  );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsBDisc")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedCvsBDisc").c_str(),  Jet_DeepFlavourPrunedCvsBDisc  );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedB")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedB"       ).c_str(),  Jet_DeepFlavourPrunedB         );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedBB")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedBB"      ).c_str(),  Jet_DeepFlavourPrunedBB        );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedLEPB")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedLEPB"    ).c_str(),  Jet_DeepFlavourPrunedLEPB      );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedC")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedC"       ).c_str(),  Jet_DeepFlavourPrunedC         );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedUDS")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedUDS"     ).c_str(),  Jet_DeepFlavourPrunedUDS       );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedG")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedG"       ).c_str(),  Jet_DeepFlavourPrunedG         );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedBDiscN")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedBDiscN"   ).c_str(), Jet_DeepFlavourPrunedBDiscN   );
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsLDiscN")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedCvsLDiscN").c_str(), Jet_DeepFlavourPrunedCvsLDiscN);
      if(variableParser.isToBeStored(name+"Jet_DeepFlavourPrunedCvsBDiscN")) tree->SetBranchAddress((name+"Jet_DeepFlavourPrunedCvsBDiscN").c_str(), Jet_DeepFlavourPrunedCvsBDiscN);

      if(variableParser.isToBeStored(name+"Jet_DeepCSVBDisc")) tree->SetBranchAddress((name+"Jet_DeepCSVBDisc"      ).c_str(), Jet_DeepCSVBDisc          );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVBDiscN")) tree->SetBranchAddress((name+"Jet_DeepCSVBDiscN"     ).c_str(), Jet_DeepCSVBDiscN        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVBDiscP")) tree->SetBranchAddress((name+"Jet_DeepCSVBDiscP"     ).c_str(), Jet_DeepCSVBDiscP        );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsLDisc")) tree->SetBranchAddress((name+"Jet_DeepCSVCvsLDisc"   ).c_str(), Jet_DeepCSVCvsLDisc    );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsLDiscN")) tree->SetBranchAddress((name+"Jet_DeepCSVCvsLDiscN"  ).c_str(), Jet_DeepCSVCvsLDiscN  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsLDiscP")) tree->SetBranchAddress((name+"Jet_DeepCSVCvsLDiscP"  ).c_str(), Jet_DeepCSVCvsLDiscP  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsBDisc")) tree->SetBranchAddress((name+"Jet_DeepCSVCvsBDisc"   ).c_str(), Jet_DeepCSVCvsBDisc    );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsBDiscN")) tree->SetBranchAddress((name+"Jet_DeepCSVCvsBDiscN"  ).c_str(), Jet_DeepCSVCvsBDiscN  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVCvsBDiscP")) tree->SetBranchAddress((name+"Jet_DeepCSVCvsBDiscP"  ).c_str(), Jet_DeepCSVCvsBDiscP  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVb")) tree->SetBranchAddress((name+"Jet_DeepCSVb"  ).c_str(),         Jet_DeepCSVb                  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVc")) tree->SetBranchAddress((name+"Jet_DeepCSVc"  ).c_str(),         Jet_DeepCSVc                  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVl")) tree->SetBranchAddress((name+"Jet_DeepCSVl"  ).c_str(),         Jet_DeepCSVl                  );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbb")) tree->SetBranchAddress((name+"Jet_DeepCSVbb" ).c_str(),         Jet_DeepCSVbb                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVcc")) tree->SetBranchAddress((name+"Jet_DeepCSVcc" ).c_str(),         Jet_DeepCSVcc                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbN")) tree->SetBranchAddress((name+"Jet_DeepCSVbN" ).c_str(),         Jet_DeepCSVbN                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVcN")) tree->SetBranchAddress((name+"Jet_DeepCSVcN" ).c_str(),         Jet_DeepCSVcN                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVlN")) tree->SetBranchAddress((name+"Jet_DeepCSVlN" ).c_str(),         Jet_DeepCSVlN                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbbN")) tree->SetBranchAddress((name+"Jet_DeepCSVbbN").c_str(),         Jet_DeepCSVbbN              );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVccN")) tree->SetBranchAddress((name+"Jet_DeepCSVccN").c_str(),         Jet_DeepCSVccN              );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbP")) tree->SetBranchAddress((name+"Jet_DeepCSVbP" ).c_str(),         Jet_DeepCSVbP                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVcP")) tree->SetBranchAddress((name+"Jet_DeepCSVcP" ).c_str(),         Jet_DeepCSVcP                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVlP")) tree->SetBranchAddress((name+"Jet_DeepCSVlP" ).c_str(),         Jet_DeepCSVlP                );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVbbP")) tree->SetBranchAddress((name+"Jet_DeepCSVbbP").c_str(),         Jet_DeepCSVbbP              );
      if(variableParser.isToBeStored(name+"Jet_DeepCSVccP")) tree->SetBranchAddress((name+"Jet_DeepCSVccP").c_str(),         Jet_DeepCSVccP              );

      //--------------------------------------
      // pf electron information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstSE")) tree->SetBranchAddress((name+"Jet_nFirstSE").c_str(),    Jet_nFirstSE     );
      if(variableParser.isToBeStored(name+"Jet_nLastSE")) tree->SetBranchAddress((name+"Jet_nLastSE").c_str(),     Jet_nLastSE       );

      if(variableParser.isToBeStored(name+"nPFElectron")) tree->SetBranchAddress((name+"nPFElectron").c_str()         ,&nPFElectron        );
      if(variableParser.isToBeStored(name+"PFElectron_IdxJet")) tree->SetBranchAddress((name+"PFElectron_IdxJet").c_str()   ,PFElectron_IdxJet   );
      if(variableParser.isToBeStored(name+"PFElectron_pt")) tree->SetBranchAddress((name+"PFElectron_pt").c_str()       ,PFElectron_pt       );
      if(variableParser.isToBeStored(name+"PFElectron_eta")) tree->SetBranchAddress((name+"PFElectron_eta").c_str()      ,PFElectron_eta      );
      if(variableParser.isToBeStored(name+"PFElectron_phi")) tree->SetBranchAddress((name+"PFElectron_phi").c_str()      ,PFElectron_phi      );
      if(variableParser.isToBeStored(name+"PFElectron_ptrel")) tree->SetBranchAddress((name+"PFElectron_ptrel").c_str()    ,PFElectron_ptrel    );
      if(variableParser.isToBeStored(name+"PFElectron_deltaR")) tree->SetBranchAddress((name+"PFElectron_deltaR").c_str()   ,PFElectron_deltaR   );
      if(variableParser.isToBeStored(name+"PFElectron_ratio")) tree->SetBranchAddress((name+"PFElectron_ratio").c_str()    ,PFElectron_ratio    );
      if(variableParser.isToBeStored(name+"PFElectron_ratioRel")) tree->SetBranchAddress((name+"PFElectron_ratioRel").c_str() ,PFElectron_ratioRel );
      if(variableParser.isToBeStored(name+"PFElectron_IP")) tree->SetBranchAddress((name+"PFElectron_IP").c_str()       ,PFElectron_IP       );
      if(variableParser.isToBeStored(name+"PFElectron_IP2D")) tree->SetBranchAddress((name+"PFElectron_IP2D").c_str()     ,PFElectron_IP2D     );

      //--------------------------------------
      // pf muon information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstSM")) tree->SetBranchAddress((name+"Jet_nFirstSM").c_str()        ,Jet_nFirstSM         );
      if(variableParser.isToBeStored(name+"Jet_nLastSM")) tree->SetBranchAddress((name+"Jet_nLastSM").c_str()         ,Jet_nLastSM           );

      if(variableParser.isToBeStored(name+"nPFMuon")) tree->SetBranchAddress((name+"nPFMuon").c_str()            ,&nPFMuon            );
      if(variableParser.isToBeStored(name+"PFMuon_IdxJet")) tree->SetBranchAddress((name+"PFMuon_IdxJet").c_str()      ,PFMuon_IdxJet       );
      if(variableParser.isToBeStored(name+"PFMuon_nMuHit")) tree->SetBranchAddress((name+"PFMuon_nMuHit").c_str()      ,PFMuon_nMuHit       );
      if(variableParser.isToBeStored(name+"PFMuon_nTkHit")) tree->SetBranchAddress((name+"PFMuon_nTkHit").c_str()      ,PFMuon_nTkHit       );
      if(variableParser.isToBeStored(name+"PFMuon_nPixHit")) tree->SetBranchAddress((name+"PFMuon_nPixHit").c_str()     ,PFMuon_nPixHit      );
      if(variableParser.isToBeStored(name+"PFMuon_nOutHit")) tree->SetBranchAddress((name+"PFMuon_nOutHit").c_str()     ,PFMuon_nOutHit      );
      if(variableParser.isToBeStored(name+"PFMuon_nTkLwM")) tree->SetBranchAddress((name+"PFMuon_nTkLwM").c_str()      ,PFMuon_nTkLwM       );
      if(variableParser.isToBeStored(name+"PFMuon_nPixLwM")) tree->SetBranchAddress((name+"PFMuon_nPixLwM").c_str()     ,PFMuon_nPixLwM      );
      if(variableParser.isToBeStored(name+"PFMuon_nMatched")) tree->SetBranchAddress((name+"PFMuon_nMatched").c_str()    ,PFMuon_nMatched     );
      if(variableParser.isToBeStored(name+"PFMuon_chi2")) tree->SetBranchAddress((name+"PFMuon_chi2").c_str()        ,PFMuon_chi2         );
      if(variableParser.isToBeStored(name+"PFMuon_chi2Tk")) tree->SetBranchAddress((name+"PFMuon_chi2Tk").c_str()      ,PFMuon_chi2Tk       );
      if(variableParser.isToBeStored(name+"PFMuon_isGlobal")) tree->SetBranchAddress((name+"PFMuon_isGlobal").c_str()    ,PFMuon_isGlobal     );
      if(variableParser.isToBeStored(name+"PFMuon_hist")) tree->SetBranchAddress((name+"PFMuon_hist").c_str()        ,PFMuon_hist         );
      if(variableParser.isToBeStored(name+"PFMuon_pt")) tree->SetBranchAddress((name+"PFMuon_pt").c_str()          ,PFMuon_pt           );
      if(variableParser.isToBeStored(name+"PFMuon_eta")) tree->SetBranchAddress((name+"PFMuon_eta").c_str()         ,PFMuon_eta          );
      if(variableParser.isToBeStored(name+"PFMuon_phi")) tree->SetBranchAddress((name+"PFMuon_phi").c_str()         ,PFMuon_phi          );
      if(variableParser.isToBeStored(name+"PFMuon_ptrel")) tree->SetBranchAddress((name+"PFMuon_ptrel").c_str()       ,PFMuon_ptrel        );
      if(variableParser.isToBeStored(name+"PFMuon_deltaR")) tree->SetBranchAddress((name+"PFMuon_deltaR").c_str()      ,PFMuon_deltaR       );
      if(variableParser.isToBeStored(name+"PFMuon_ratio")) tree->SetBranchAddress((name+"PFMuon_ratio").c_str()       ,PFMuon_ratio        );
      if(variableParser.isToBeStored(name+"PFMuon_ratioRel")) tree->SetBranchAddress((name+"PFMuon_ratioRel").c_str()    ,PFMuon_ratioRel     );
      if(variableParser.isToBeStored(name+"PFMuon_IP")) tree->SetBranchAddress((name+"PFMuon_IP").c_str()          ,PFMuon_IP           );
      if(variableParser.isToBeStored(name+"PFMuon_IP2D")) tree->SetBranchAddress((name+"PFMuon_IP2D").c_str()        ,PFMuon_IP2D         );
      if(variableParser.isToBeStored(name+"PFMuon_IPsig")) tree->SetBranchAddress((name+"PFMuon_IPsig").c_str()       ,PFMuon_IPsig        );
      if(variableParser.isToBeStored(name+"PFMuon_IP2Dsig")) tree->SetBranchAddress((name+"PFMuon_IP2Dsig").c_str()     ,PFMuon_IP2Dsig      );
      if(variableParser.isToBeStored(name+"PFMuon_dz")) tree->SetBranchAddress((name+"PFMuon_dz").c_str()          ,PFMuon_dz           );
      if(variableParser.isToBeStored(name+"PFMuon_GoodQuality")) tree->SetBranchAddress((name+"PFMuon_GoodQuality").c_str() ,PFMuon_GoodQuality  );

      //--------------------------------------
      // secondary vertex information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstSV")) tree->SetBranchAddress((name+"Jet_nFirstSV").c_str(),       Jet_nFirstSV    );
      if(variableParser.isToBeStored(name+"Jet_nLastSV")) tree->SetBranchAddress((name+"Jet_nLastSV").c_str(),        Jet_nLastSV     );

      if(variableParser.isToBeStored(name+"nSV")) tree->SetBranchAddress((name+"nSV").c_str()                ,&nSV               );
      if(variableParser.isToBeStored(name+"SV_x")) tree->SetBranchAddress((name+"SV_x").c_str()               ,SV_x               );
      if(variableParser.isToBeStored(name+"SV_y")) tree->SetBranchAddress((name+"SV_y").c_str()               ,SV_y               );
      if(variableParser.isToBeStored(name+"SV_z")) tree->SetBranchAddress((name+"SV_z").c_str()               ,SV_z               );
      if(variableParser.isToBeStored(name+"SV_ex")) tree->SetBranchAddress((name+"SV_ex").c_str()              ,SV_ex              );
      if(variableParser.isToBeStored(name+"SV_ey")) tree->SetBranchAddress((name+"SV_ey").c_str()              ,SV_ey              );
      if(variableParser.isToBeStored(name+"SV_ez")) tree->SetBranchAddress((name+"SV_ez").c_str()              ,SV_ez              );
      if(variableParser.isToBeStored(name+"SV_chi2")) tree->SetBranchAddress((name+"SV_chi2").c_str()            ,SV_chi2            );
      if(variableParser.isToBeStored(name+"SV_ndf")) tree->SetBranchAddress((name+"SV_ndf").c_str()             ,SV_ndf             );
      if(variableParser.isToBeStored(name+"SV_flight")) tree->SetBranchAddress((name+"SV_flight").c_str()          ,SV_flight          );
      if(variableParser.isToBeStored(name+"SV_flightErr")) tree->SetBranchAddress((name+"SV_flightErr").c_str()       ,SV_flightErr       );
      if(variableParser.isToBeStored(name+"SV_deltaR_jet")) tree->SetBranchAddress((name+"SV_deltaR_jet").c_str()      ,SV_deltaR_jet      );
      if(variableParser.isToBeStored(name+"SV_deltaR_sum_jet")) tree->SetBranchAddress((name+"SV_deltaR_sum_jet").c_str()  ,SV_deltaR_sum_jet  );
      if(variableParser.isToBeStored(name+"SV_deltaR_sum_dir")) tree->SetBranchAddress((name+"SV_deltaR_sum_dir").c_str()  ,SV_deltaR_sum_dir  );
      if(variableParser.isToBeStored(name+"SV_vtx_pt")) tree->SetBranchAddress((name+"SV_vtx_pt").c_str()          ,SV_vtx_pt          );
      if(variableParser.isToBeStored(name+"SV_flight2D")) tree->SetBranchAddress((name+"SV_flight2D").c_str()        ,SV_flight2D        );
      if(variableParser.isToBeStored(name+"SV_flight2DErr")) tree->SetBranchAddress((name+"SV_flight2DErr").c_str()     ,SV_flight2DErr     );
      if(variableParser.isToBeStored(name+"SV_totCharge")) tree->SetBranchAddress((name+"SV_totCharge").c_str()       ,SV_totCharge       );
      if(variableParser.isToBeStored(name+"SV_vtxDistJetAxis")) tree->SetBranchAddress((name+"SV_vtxDistJetAxis").c_str()  ,SV_vtxDistJetAxis  );
      if(variableParser.isToBeStored(name+"SV_nTrk")) tree->SetBranchAddress((name+"SV_nTrk").c_str()            ,SV_nTrk            );
      if(variableParser.isToBeStored(name+"SV_mass")) tree->SetBranchAddress((name+"SV_mass").c_str()            ,SV_mass            );
      if(variableParser.isToBeStored(name+"SV_vtx_eta")) tree->SetBranchAddress((name+"SV_vtx_eta").c_str()         ,SV_vtx_eta         );
      if(variableParser.isToBeStored(name+"SV_vtx_phi")) tree->SetBranchAddress((name+"SV_vtx_phi").c_str()         ,SV_vtx_phi         );
      if(variableParser.isToBeStored(name+"SV_EnergyRatio")) tree->SetBranchAddress((name+"SV_EnergyRatio").c_str()     ,SV_EnergyRatio     );
      if(variableParser.isToBeStored(name+"SV_dir_x")) tree->SetBranchAddress((name+"SV_dir_x").c_str()           ,SV_dir_x           );
      if(variableParser.isToBeStored(name+"SV_dir_y")) tree->SetBranchAddress((name+"SV_dir_y").c_str()           ,SV_dir_y           );
      if(variableParser.isToBeStored(name+"SV_dir_z")) tree->SetBranchAddress((name+"SV_dir_z").c_str()           ,SV_dir_z           );

      //--------------------------------------
      // track information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrack")) tree->SetBranchAddress((name+"Jet_nFirstTrack").c_str(),  Jet_nFirstTrack );
      if(variableParser.isToBeStored(name+"Jet_nLastTrack")) tree->SetBranchAddress((name+"Jet_nLastTrack").c_str(),   Jet_nLastTrack  );
      if(variableParser.isToBeStored(name+"nTrack")) tree->SetBranchAddress((name+"nTrack").c_str()           ,&nTrack          );
      if(variableParser.isToBeStored(name+"Track_dxy")) tree->SetBranchAddress((name+"Track_dxy").c_str()        ,Track_dxy             );
      if(variableParser.isToBeStored(name+"Track_dz")) tree->SetBranchAddress((name+"Track_dz").c_str()         ,Track_dz         );
      if(variableParser.isToBeStored(name+"Track_dxyError")) tree->SetBranchAddress((name+"Track_dxyError").c_str()   ,Track_dxyError   );
      if(variableParser.isToBeStored(name+"Track_dzError")) tree->SetBranchAddress((name+"Track_dzError").c_str()    ,Track_dzError    );
      if(variableParser.isToBeStored(name+"Track_sign2D")) tree->SetBranchAddress((name+"Track_sign2D").c_str()     ,Track_sign2D    );
      if(variableParser.isToBeStored(name+"Track_sign3D")) tree->SetBranchAddress((name+"Track_sign3D").c_str()     ,Track_sign3D    );
      if(variableParser.isToBeStored(name+"Track_length")) tree->SetBranchAddress((name+"Track_length").c_str()     ,Track_length     );
      if(variableParser.isToBeStored(name+"Track_dist")) tree->SetBranchAddress((name+"Track_dist").c_str()       ,Track_dist            );
      if(variableParser.isToBeStored(name+"Track_IP2D")) tree->SetBranchAddress((name+"Track_IP2D").c_str()       ,Track_IP2D            );
      if(variableParser.isToBeStored(name+"Track_IP2Dsig")) tree->SetBranchAddress((name+"Track_IP2Dsig").c_str()    ,Track_IP2Dsig    );
      if(variableParser.isToBeStored(name+"Track_IP2Derr")) tree->SetBranchAddress((name+"Track_IP2Derr").c_str()    ,Track_IP2Derr    );
      if(variableParser.isToBeStored(name+"Track_IP")) tree->SetBranchAddress((name+"Track_IP").c_str()         ,Track_IP         );
      if(variableParser.isToBeStored(name+"Track_IPsig")) tree->SetBranchAddress((name+"Track_IPsig").c_str()      ,Track_IPsig      );
      if(variableParser.isToBeStored(name+"Track_IPerr")) tree->SetBranchAddress((name+"Track_IPerr").c_str()      ,Track_IPerr      );
      if(variableParser.isToBeStored(name+"Track_Proba")) tree->SetBranchAddress((name+"Track_Proba").c_str()      ,Track_Proba      );
      if(variableParser.isToBeStored(name+"Track_p")) tree->SetBranchAddress((name+"Track_p").c_str()          ,Track_p          );
      if(variableParser.isToBeStored(name+"Track_pt")) tree->SetBranchAddress((name+"Track_pt").c_str()         ,Track_pt         );
      if(variableParser.isToBeStored(name+"Track_eta")) tree->SetBranchAddress((name+"Track_eta").c_str()        ,Track_eta             );
      if(variableParser.isToBeStored(name+"Track_phi")) tree->SetBranchAddress((name+"Track_phi").c_str()        ,Track_phi             );
      if(variableParser.isToBeStored(name+"Track_chi2")) tree->SetBranchAddress((name+"Track_chi2").c_str()       ,Track_chi2            );
      if(variableParser.isToBeStored(name+"Track_charge")) tree->SetBranchAddress((name+"Track_charge").c_str()     ,Track_charge     );
      if(variableParser.isToBeStored(name+"Track_history")) tree->SetBranchAddress((name+"Track_history").c_str()    ,Track_history    );
      if(variableParser.isToBeStored(name+"Track_nHitStrip")) tree->SetBranchAddress((name+"Track_nHitStrip").c_str()  ,Track_nHitStrip  );
      if(variableParser.isToBeStored(name+"Track_nHitPixel")) tree->SetBranchAddress((name+"Track_nHitPixel").c_str()  ,Track_nHitPixel  );
      if(variableParser.isToBeStored(name+"Track_nHitAll")) tree->SetBranchAddress((name+"Track_nHitAll").c_str()    ,Track_nHitAll    );
      if(variableParser.isToBeStored(name+"Track_nHitTIB")) tree->SetBranchAddress((name+"Track_nHitTIB").c_str()    ,Track_nHitTIB    );
      if(variableParser.isToBeStored(name+"Track_nHitTID")) tree->SetBranchAddress((name+"Track_nHitTID").c_str()    ,Track_nHitTID    );
      if(variableParser.isToBeStored(name+"Track_nHitTOB")) tree->SetBranchAddress((name+"Track_nHitTOB").c_str()    ,Track_nHitTOB    );
      if(variableParser.isToBeStored(name+"Track_nHitTEC")) tree->SetBranchAddress((name+"Track_nHitTEC").c_str()    ,Track_nHitTEC    );
      if(variableParser.isToBeStored(name+"Track_nHitPXB")) tree->SetBranchAddress((name+"Track_nHitPXB").c_str()    ,Track_nHitPXB    );
      if(variableParser.isToBeStored(name+"Track_nHitPXF")) tree->SetBranchAddress((name+"Track_nHitPXF").c_str()    ,Track_nHitPXF    );
      if(variableParser.isToBeStored(name+"Track_isHitL1")) tree->SetBranchAddress((name+"Track_isHitL1").c_str()    ,Track_isHitL1    );
      if(variableParser.isToBeStored(name+"Track_PV")) tree->SetBranchAddress((name+"Track_PV").c_str()         ,Track_PV         );
      if(variableParser.isToBeStored(name+"Track_SV")) tree->SetBranchAddress((name+"Track_SV").c_str()         ,Track_SV         );
      if(variableParser.isToBeStored(name+"Track_PVweight")) tree->SetBranchAddress((name+"Track_PVweight").c_str()   ,Track_PVweight   );
      if(variableParser.isToBeStored(name+"Track_SVweight")) tree->SetBranchAddress((name+"Track_SVweight").c_str()   ,Track_SVweight   );
      if(variableParser.isToBeStored(name+"Track_isfromSV")) tree->SetBranchAddress((name+"Track_isfromSV").c_str()   ,Track_isfromSV   );
      if(variableParser.isToBeStored(name+"Track_isfromV0")) tree->SetBranchAddress((name+"Track_isfromV0").c_str()   ,Track_isfromV0   );
      if(variableParser.isToBeStored(name+"Track_category")) tree->SetBranchAddress((name+"Track_category").c_str()   ,Track_category   );

      //--------------------------------------
      // track truth information
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Track_TPAssociationQuality")) tree->SetBranchAddress((name+"Track_TPAssociationQuality").c_str()          ,Track_TPAssociationQuality             );
      if(variableParser.isToBeStored(name+"Track_idxMatchedTP")) tree->SetBranchAddress((name+"Track_idxMatchedTP").c_str()                  ,Track_idxMatchedTP                     );

      if(variableParser.isToBeStored(name+"Jet_nFirstTrackTruth")) tree->SetBranchAddress((name+"Jet_nFirstTrackTruth").c_str()                ,Jet_nFirstTrackTruth                   );
      if(variableParser.isToBeStored(name+"Jet_nLastTrackTruth")) tree->SetBranchAddress((name+"Jet_nLastTrackTruth").c_str()                 ,Jet_nLastTrackTruth                    );

      if(variableParser.isToBeStored(name+"nTrackTruth")) tree->SetBranchAddress((name+"nTrackTruth").c_str()                ,&nTrackTruth                           );
      if(variableParser.isToBeStored(name+"TrackTruth_idxMatchedTrack")) tree->SetBranchAddress((name+"TrackTruth_idxMatchedTrack").c_str()          ,TrackTruth_idxMatchedTrack             );
      if(variableParser.isToBeStored(name+"TrackTruth_p")) tree->SetBranchAddress((name+"TrackTruth_p").c_str()                        ,TrackTruth_p                           );
      if(variableParser.isToBeStored(name+"TrackTruth_pt")) tree->SetBranchAddress((name+"TrackTruth_pt").c_str()                       ,TrackTruth_pt                          );
      if(variableParser.isToBeStored(name+"TrackTruth_eta")) tree->SetBranchAddress((name+"TrackTruth_eta").c_str()                      ,TrackTruth_eta                         );
      if(variableParser.isToBeStored(name+"TrackTruth_phi")) tree->SetBranchAddress((name+"TrackTruth_phi").c_str()                      ,TrackTruth_phi                         );
      if(variableParser.isToBeStored(name+"TrackTruth_charge")) tree->SetBranchAddress((name+"TrackTruth_charge").c_str()                   ,TrackTruth_charge                      );
      if(variableParser.isToBeStored(name+"TrackTruth_pdgid")) tree->SetBranchAddress((name+"TrackTruth_pdgid").c_str()                    ,TrackTruth_pdgid                       );
      if(variableParser.isToBeStored(name+"TrackTruth_dxy")) tree->SetBranchAddress((name+"TrackTruth_dxy").c_str()                      ,TrackTruth_dxy                         );
      if(variableParser.isToBeStored(name+"TrackTruth_dz")) tree->SetBranchAddress((name+"TrackTruth_dz").c_str()                       ,TrackTruth_dz                          );
      if(variableParser.isToBeStored(name+"TrackTruth_nHitAll")) tree->SetBranchAddress((name+"TrackTruth_nHitAll").c_str()                  ,TrackTruth_nHitAll                     );
      if(variableParser.isToBeStored(name+"TrackTruth_nHitPixel")) tree->SetBranchAddress((name+"TrackTruth_nHitPixel").c_str()                ,TrackTruth_nHitPixel                   );
      if(variableParser.isToBeStored(name+"TrackTruth_nHitStrip")) tree->SetBranchAddress((name+"TrackTruth_nHitStrip").c_str()                ,TrackTruth_nHitStrip                   );

      //--------------------------------------
      // Inclusive Track information for PtRel template
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkInc")) tree->SetBranchAddress((name+"Jet_nFirstTrkInc").c_str(), Jet_nFirstTrkInc );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkInc")) tree->SetBranchAddress((name+"Jet_nLastTrkInc").c_str(),  Jet_nLastTrkInc  );

      if(variableParser.isToBeStored(name+"nTrkInc")) tree->SetBranchAddress((name+"nTrkInc").c_str()      ,&nTrkInc     );
      if(variableParser.isToBeStored(name+"TrkInc_pt")) tree->SetBranchAddress((name+"TrkInc_pt").c_str()    ,TrkInc_pt    );
      if(variableParser.isToBeStored(name+"TrkInc_eta")) tree->SetBranchAddress((name+"TrkInc_eta").c_str()   ,TrkInc_eta   );
      if(variableParser.isToBeStored(name+"TrkInc_phi")) tree->SetBranchAddress((name+"TrkInc_phi").c_str()   ,TrkInc_phi   );
      if(variableParser.isToBeStored(name+"TrkInc_ptrel")) tree->SetBranchAddress((name+"TrkInc_ptrel").c_str() ,TrkInc_ptrel );
      if(variableParser.isToBeStored(name+"TrkInc_IPsig")) tree->SetBranchAddress((name+"TrkInc_IPsig").c_str() ,TrkInc_IPsig );
      if(variableParser.isToBeStored(name+"TrkInc_IP")) tree->SetBranchAddress((name+"TrkInc_IP").c_str()    ,TrkInc_IP    );

      //--------------------------------------
      // TagInfo TaggingVariables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkTagVar")) tree->SetBranchAddress((name+"Jet_nFirstTrkTagVar").c_str()  ,Jet_nFirstTrkTagVar   );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkTagVar")) tree->SetBranchAddress((name+"Jet_nLastTrkTagVar").c_str()   ,Jet_nLastTrkTagVar         );
      if(variableParser.isToBeStored(name+"Jet_nFirstSVTagVar")) tree->SetBranchAddress((name+"Jet_nFirstSVTagVar").c_str()   ,Jet_nFirstSVTagVar     );
      if(variableParser.isToBeStored(name+"Jet_nLastSVTagVar")) tree->SetBranchAddress((name+"Jet_nLastSVTagVar").c_str()    ,Jet_nLastSVTagVar       );

      if(variableParser.isToBeStored(name+"TagVar_jetNTracks")) tree->SetBranchAddress((name+"TagVar_jetNTracks").c_str()                   ,TagVar_jetNTracks                                      );
      if(variableParser.isToBeStored(name+"TagVar_jetNSecondaryVertices")) tree->SetBranchAddress((name+"TagVar_jetNSecondaryVertices").c_str()        ,TagVar_jetNSecondaryVertices                );
      //run 2 id: https://twiki.cern.ch/CMS/JetID13TeVRun2017
      if(variableParser.isToBeStored(name+"TagVar_chargedHadronEnergyFraction")) tree->SetBranchAddress((name+"TagVar_chargedHadronEnergyFraction").c_str()  ,TagVar_chargedHadronEnergyFraction    );
      if(variableParser.isToBeStored(name+"TagVar_neutralHadronEnergyFraction")) tree->SetBranchAddress((name+"TagVar_neutralHadronEnergyFraction").c_str()  ,TagVar_neutralHadronEnergyFraction    );
      if(variableParser.isToBeStored(name+"TagVar_photonEnergyFraction")) tree->SetBranchAddress((name+"TagVar_photonEnergyFraction").c_str()         ,TagVar_photonEnergyFraction                  );
      if(variableParser.isToBeStored(name+"TagVar_electronEnergyFraction")) tree->SetBranchAddress((name+"TagVar_electronEnergyFraction").c_str()       ,TagVar_electronEnergyFraction              );
      if(variableParser.isToBeStored(name+"TagVar_muonEnergyFraction")) tree->SetBranchAddress((name+"TagVar_muonEnergyFraction").c_str()           ,TagVar_muonEnergyFraction                      );
      if(variableParser.isToBeStored(name+"TagVar_chargedHadronMultiplicity")) tree->SetBranchAddress((name+"TagVar_chargedHadronMultiplicity").c_str()    ,TagVar_chargedHadronMultiplicity        );
      if(variableParser.isToBeStored(name+"TagVar_neutralHadronMultiplicity")) tree->SetBranchAddress((name+"TagVar_neutralHadronMultiplicity").c_str()    ,TagVar_neutralHadronMultiplicity        );
      if(variableParser.isToBeStored(name+"TagVar_photonMultiplicity")) tree->SetBranchAddress((name+"TagVar_photonMultiplicity").c_str()           ,TagVar_photonMultiplicity                      );
      if(variableParser.isToBeStored(name+"TagVar_electronMultiplicity")) tree->SetBranchAddress((name+"TagVar_electronMultiplicity").c_str()         ,TagVar_electronMultiplicity                  );
      if(variableParser.isToBeStored(name+"TagVar_muonMultiplicity")) tree->SetBranchAddress((name+"TagVar_muonMultiplicity").c_str()             ,TagVar_muonMultiplicity                          );
      if(variableParser.isToBeStored(name+"TagVar_neutralEmEnergyFraction")) tree->SetBranchAddress((name+"TagVar_neutralEmEnergyFraction").c_str()      ,TagVar_neutralEmEnergyFraction            );
      if(variableParser.isToBeStored(name+"TagVar_numberOfDaughters")) tree->SetBranchAddress((name+"TagVar_numberOfDaughters").c_str()            ,TagVar_numberOfDaughters                        );
      if(variableParser.isToBeStored(name+"TagVar_chargedMultiplicity")) tree->SetBranchAddress((name+"TagVar_chargedMultiplicity").c_str()          ,TagVar_chargedMultiplicity                    );
      if(variableParser.isToBeStored(name+"TagVar_chargedEmEnergyFraction")) tree->SetBranchAddress((name+"TagVar_chargedEmEnergyFraction").c_str()      ,TagVar_chargedEmEnergyFraction            );
      if(variableParser.isToBeStored(name+"TagVar_neutralMultiplicity")) tree->SetBranchAddress((name+"TagVar_neutralMultiplicity").c_str()          ,TagVar_neutralMultiplicity                    );

      if(variableParser.isToBeStored(name+"nTrkTagVar")) tree->SetBranchAddress((name+"nTrkTagVar").c_str()               ,&nTrkTagVar                                                );
      if(variableParser.isToBeStored(name+"TagVar_trackMomentum")) tree->SetBranchAddress((name+"TagVar_trackMomentum").c_str()     ,TagVar_trackMomentum         );
      if(variableParser.isToBeStored(name+"TagVar_trackEta")) tree->SetBranchAddress((name+"TagVar_trackEta").c_str()          ,TagVar_trackEta                   );
      if(variableParser.isToBeStored(name+"TagVar_trackPhi")) tree->SetBranchAddress((name+"TagVar_trackPhi").c_str()          ,TagVar_trackPhi                   );
      if(variableParser.isToBeStored(name+"TagVar_trackPtRel")) tree->SetBranchAddress((name+"TagVar_trackPtRel").c_str()        ,TagVar_trackPtRel               );
      if(variableParser.isToBeStored(name+"TagVar_trackPPar")) tree->SetBranchAddress((name+"TagVar_trackPPar").c_str()         ,TagVar_trackPPar                 );
      if(variableParser.isToBeStored(name+"TagVar_trackEtaRel")) tree->SetBranchAddress((name+"TagVar_trackEtaRel").c_str()       ,TagVar_trackEtaRel             );
      if(variableParser.isToBeStored(name+"TagVar_trackDeltaR")) tree->SetBranchAddress((name+"TagVar_trackDeltaR").c_str()       ,TagVar_trackDeltaR             );
      if(variableParser.isToBeStored(name+"TagVar_trackPtRatio")) tree->SetBranchAddress((name+"TagVar_trackPtRatio").c_str()      ,TagVar_trackPtRatio           );
      if(variableParser.isToBeStored(name+"TagVar_trackPParRatio")) tree->SetBranchAddress((name+"TagVar_trackPParRatio").c_str()    ,TagVar_trackPParRatio       );
      if(variableParser.isToBeStored(name+"TagVar_trackSip2dVal")) tree->SetBranchAddress((name+"TagVar_trackSip2dVal").c_str()     ,TagVar_trackSip2dVal         );
      if(variableParser.isToBeStored(name+"TagVar_trackSip2dSig")) tree->SetBranchAddress((name+"TagVar_trackSip2dSig").c_str()     ,TagVar_trackSip2dSig         );
      if(variableParser.isToBeStored(name+"TagVar_trackSip3dVal")) tree->SetBranchAddress((name+"TagVar_trackSip3dVal").c_str()     ,TagVar_trackSip3dVal         );
      if(variableParser.isToBeStored(name+"TagVar_trackSip3dSig")) tree->SetBranchAddress((name+"TagVar_trackSip3dSig").c_str()     ,TagVar_trackSip3dSig         );
      if(variableParser.isToBeStored(name+"TagVar_trackDecayLenVal")) tree->SetBranchAddress((name+"TagVar_trackDecayLenVal").c_str()  ,TagVar_trackDecayLenVal   );
      if(variableParser.isToBeStored(name+"TagVar_trackDecayLenSig")) tree->SetBranchAddress((name+"TagVar_trackDecayLenSig").c_str()  ,TagVar_trackDecayLenSig   );
      if(variableParser.isToBeStored(name+"TagVar_trackJetDistVal")) tree->SetBranchAddress((name+"TagVar_trackJetDistVal").c_str()   ,TagVar_trackJetDistVal     );
      if(variableParser.isToBeStored(name+"TagVar_trackJetDistSig")) tree->SetBranchAddress((name+"TagVar_trackJetDistSig").c_str()   ,TagVar_trackJetDistSig     );
      if(variableParser.isToBeStored(name+"TagVar_trackChi2")) tree->SetBranchAddress((name+"TagVar_trackChi2").c_str()         ,TagVar_trackChi2                 );
      if(variableParser.isToBeStored(name+"TagVar_trackNTotalHits")) tree->SetBranchAddress((name+"TagVar_trackNTotalHits").c_str()   ,TagVar_trackNTotalHits     );
      if(variableParser.isToBeStored(name+"TagVar_trackNPixelHits")) tree->SetBranchAddress((name+"TagVar_trackNPixelHits").c_str()   ,TagVar_trackNPixelHits     );

      if(variableParser.isToBeStored(name+"nSVTagVar")) tree->SetBranchAddress((name+"nSVTagVar").c_str()                       ,&nSVTagVar                                                           );
      if(variableParser.isToBeStored(name+"TagVar_vertexMass")) tree->SetBranchAddress((name+"TagVar_vertexMass").c_str()               ,TagVar_vertexMass                         );
      if(variableParser.isToBeStored(name+"TagVar_vertexNTracks")) tree->SetBranchAddress((name+"TagVar_vertexNTracks").c_str()            ,TagVar_vertexNTracks                   );
      if(variableParser.isToBeStored(name+"TagVar_vertexJetDeltaR")) tree->SetBranchAddress((name+"TagVar_vertexJetDeltaR").c_str()          ,TagVar_vertexJetDeltaR               );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance2dVal")) tree->SetBranchAddress((name+"TagVar_flightDistance2dVal").c_str()      ,TagVar_flightDistance2dVal       );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance2dSig")) tree->SetBranchAddress((name+"TagVar_flightDistance2dSig").c_str()      ,TagVar_flightDistance2dSig       );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance3dVal")) tree->SetBranchAddress((name+"TagVar_flightDistance3dVal").c_str()      ,TagVar_flightDistance3dVal       );
      if(variableParser.isToBeStored(name+"TagVar_flightDistance3dSig")) tree->SetBranchAddress((name+"TagVar_flightDistance3dSig").c_str()      ,TagVar_flightDistance3dSig       );

      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_Sip3dVal")) tree->SetBranchAddress((name+"DeepFlavourInput_charged_Sip3dVal").c_str(), DeepFlavourInput_charged_Sip3dVal);
      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_Sip3dSig")) tree->SetBranchAddress((name+"DeepFlavourInput_charged_Sip3dSig").c_str(), DeepFlavourInput_charged_Sip3dSig);
      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_quality")) tree->SetBranchAddress((name+"DeepFlavourInput_charged_quality" ).c_str(), DeepFlavourInput_charged_quality );
      if(variableParser.isToBeStored(name+"DeepFlavourInput_charged_chi2")) tree->SetBranchAddress((name+"DeepFlavourInput_charged_chi2"    ).c_str(), DeepFlavourInput_charged_chi2    );

      if(variableParser.isToBeStored(name+"DeepFlavourInput_neutral_drminsv")) tree->SetBranchAddress((name+"DeepFlavourInput_neutral_drminsv").c_str(), DeepFlavourInput_neutral_drminsv);
      if(variableParser.isToBeStored(name+"DeepFlavourInput_neutral_hadFrac")) tree->SetBranchAddress((name+"DeepFlavourInput_neutral_hadFrac").c_str(), DeepFlavourInput_neutral_hadFrac);
      if(variableParser.isToBeStored(name+"DeepFlavourInput_neutral_ptrel")) tree->SetBranchAddress((name+"DeepFlavourInput_neutral_ptrel  ").c_str(), DeepFlavourInput_neutral_ptrel  );

      if(variableParser.isToBeStored(name+"DeepFlavourInput_sv_d3d")) tree->SetBranchAddress((name+"DeepFlavourInput_sv_d3d"     ).c_str(), DeepFlavourInput_sv_d3d     );
      if(variableParser.isToBeStored(name+"DeepFlavourInput_sv_d3dsig")) tree->SetBranchAddress((name+"DeepFlavourInput_sv_d3dsig"  ).c_str(), DeepFlavourInput_sv_d3dsig  );
      if(variableParser.isToBeStored(name+"DeepFlavourInput_sv_normchi2")) tree->SetBranchAddress((name+"DeepFlavourInput_sv_normchi2").c_str(), DeepFlavourInput_sv_normchi2);

      //--------------------------------------
      // CSV TaggingVariables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"TagVarCSV_trackJetPt")) tree->SetBranchAddress((name+"TagVarCSV_trackJetPt").c_str()               ,TagVarCSV_trackJetPt                             );
      if(variableParser.isToBeStored(name+"TagVarCSV_jetNTracks")) tree->SetBranchAddress((name+"TagVarCSV_jetNTracks").c_str()               ,TagVarCSV_jetNTracks                             );
      if(variableParser.isToBeStored(name+"TagVarCSV_jetNTracksEtaRel")) tree->SetBranchAddress((name+"TagVarCSV_jetNTracksEtaRel").c_str()         ,TagVarCSV_jetNTracksEtaRel                 );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSumJetEtRatio")) tree->SetBranchAddress((name+"TagVarCSV_trackSumJetEtRatio").c_str()       ,TagVarCSV_trackSumJetEtRatio             );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSumJetDeltaR")) tree->SetBranchAddress((name+"TagVarCSV_trackSumJetDeltaR").c_str()        ,TagVarCSV_trackSumJetDeltaR               );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dValAboveCharm")) tree->SetBranchAddress((name+"TagVarCSV_trackSip2dValAboveCharm").c_str()  ,TagVarCSV_trackSip2dValAboveCharm   );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dSigAboveCharm")) tree->SetBranchAddress((name+"TagVarCSV_trackSip2dSigAboveCharm").c_str()  ,TagVarCSV_trackSip2dSigAboveCharm   );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dValAboveCharm")) tree->SetBranchAddress((name+"TagVarCSV_trackSip3dValAboveCharm").c_str()  ,TagVarCSV_trackSip3dValAboveCharm   );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dSigAboveCharm")) tree->SetBranchAddress((name+"TagVarCSV_trackSip3dSigAboveCharm").c_str()  ,TagVarCSV_trackSip3dSigAboveCharm   );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexCategory")) tree->SetBranchAddress((name+"TagVarCSV_vertexCategory").c_str()           ,TagVarCSV_vertexCategory                     );
      if(variableParser.isToBeStored(name+"TagVarCSV_jetNSecondaryVertices")) tree->SetBranchAddress((name+"TagVarCSV_jetNSecondaryVertices").c_str()    ,TagVarCSV_jetNSecondaryVertices       );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexMass")) tree->SetBranchAddress((name+"TagVarCSV_vertexMass").c_str()               ,TagVarCSV_vertexMass                             );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexNTracks")) tree->SetBranchAddress((name+"TagVarCSV_vertexNTracks").c_str()            ,TagVarCSV_vertexNTracks                       );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexEnergyRatio")) tree->SetBranchAddress((name+"TagVarCSV_vertexEnergyRatio").c_str()        ,TagVarCSV_vertexEnergyRatio               );
      if(variableParser.isToBeStored(name+"TagVarCSV_vertexJetDeltaR")) tree->SetBranchAddress((name+"TagVarCSV_vertexJetDeltaR").c_str()          ,TagVarCSV_vertexJetDeltaR                   );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance2dVal")) tree->SetBranchAddress((name+"TagVarCSV_flightDistance2dVal").c_str()      ,TagVarCSV_flightDistance2dVal           );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance2dSig")) tree->SetBranchAddress((name+"TagVarCSV_flightDistance2dSig").c_str()      ,TagVarCSV_flightDistance2dSig           );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance3dVal")) tree->SetBranchAddress((name+"TagVarCSV_flightDistance3dVal").c_str()      ,TagVarCSV_flightDistance3dVal           );
      if(variableParser.isToBeStored(name+"TagVarCSV_flightDistance3dSig")) tree->SetBranchAddress((name+"TagVarCSV_flightDistance3dSig").c_str()      ,TagVarCSV_flightDistance3dSig           );

      //--------------------------------------
      // CSV Tagging Track Variables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkTagVarCSV")) tree->SetBranchAddress((name+"Jet_nFirstTrkTagVarCSV").c_str()        ,Jet_nFirstTrkTagVarCSV               );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkTagVarCSV")) tree->SetBranchAddress((name+"Jet_nLastTrkTagVarCSV").c_str()         ,Jet_nLastTrkTagVarCSV                 );
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkEtaRelTagVarCSV")) tree->SetBranchAddress((name+"Jet_nFirstTrkEtaRelTagVarCSV").c_str()  ,Jet_nFirstTrkEtaRelTagVarCSV   );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkEtaRelTagVarCSV")) tree->SetBranchAddress((name+"Jet_nLastTrkEtaRelTagVarCSV").c_str()   ,Jet_nLastTrkEtaRelTagVarCSV         );

      if(variableParser.isToBeStored(name+"nTrkTagVarCSV")) tree->SetBranchAddress((name+"nTrkTagVarCSV").c_str()               ,&nTrkTagVarCSV                                                    );
      if(variableParser.isToBeStored(name+"nTrkEtaRelTagVarCSV")) tree->SetBranchAddress((name+"nTrkEtaRelTagVarCSV").c_str()         ,&nTrkEtaRelTagVarCSV                                        );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackMomentum")) tree->SetBranchAddress((name+"TagVarCSV_trackMomentum").c_str()     ,TagVarCSV_trackMomentum          );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackEta")) tree->SetBranchAddress((name+"TagVarCSV_trackEta").c_str()          ,TagVarCSV_trackEta                    );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPhi")) tree->SetBranchAddress((name+"TagVarCSV_trackPhi").c_str()          ,TagVarCSV_trackPhi                    );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPtRel")) tree->SetBranchAddress((name+"TagVarCSV_trackPtRel").c_str()        ,TagVarCSV_trackPtRel                );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPPar")) tree->SetBranchAddress((name+"TagVarCSV_trackPPar").c_str()         ,TagVarCSV_trackPPar                  );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackDeltaR")) tree->SetBranchAddress((name+"TagVarCSV_trackDeltaR").c_str()       ,TagVarCSV_trackDeltaR              );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPtRatio")) tree->SetBranchAddress((name+"TagVarCSV_trackPtRatio").c_str()      ,TagVarCSV_trackPtRatio            );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackPParRatio")) tree->SetBranchAddress((name+"TagVarCSV_trackPParRatio").c_str()    ,TagVarCSV_trackPParRatio        );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dVal")) tree->SetBranchAddress((name+"TagVarCSV_trackSip2dVal").c_str()     ,TagVarCSV_trackSip2dVal          );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip2dSig")) tree->SetBranchAddress((name+"TagVarCSV_trackSip2dSig").c_str()     ,TagVarCSV_trackSip2dSig          );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dVal")) tree->SetBranchAddress((name+"TagVarCSV_trackSip3dVal").c_str()     ,TagVarCSV_trackSip3dVal          );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackSip3dSig")) tree->SetBranchAddress((name+"TagVarCSV_trackSip3dSig").c_str()     ,TagVarCSV_trackSip3dSig          );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackDecayLenVal")) tree->SetBranchAddress((name+"TagVarCSV_trackDecayLenVal").c_str()  ,TagVarCSV_trackDecayLenVal    );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackDecayLenSig")) tree->SetBranchAddress((name+"TagVarCSV_trackDecayLenSig").c_str()  ,TagVarCSV_trackDecayLenSig    );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackJetDistVal")) tree->SetBranchAddress((name+"TagVarCSV_trackJetDistVal").c_str()   ,TagVarCSV_trackJetDistVal      );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackJetDistSig")) tree->SetBranchAddress((name+"TagVarCSV_trackJetDistSig").c_str()   ,TagVarCSV_trackJetDistSig      );
      if(variableParser.isToBeStored(name+"TagVarCSV_trackEtaRel")) tree->SetBranchAddress((name+"TagVarCSV_trackEtaRel").c_str()       ,TagVarCSV_trackEtaRel        );

      //--------------------------------------
      // CTag TaggingVariables
      //--------------------------------------
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsB")) tree->SetBranchAddress((name+"CTag_Jet_CvsB").c_str(),     CTag_Jet_CvsB     );
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsBN")) tree->SetBranchAddress((name+"CTag_Jet_CvsBN").c_str(),    CTag_Jet_CvsBN    );
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsBP")) tree->SetBranchAddress((name+"CTag_Jet_CvsBP").c_str(),    CTag_Jet_CvsBP    );
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsL")) tree->SetBranchAddress((name+"CTag_Jet_CvsL").c_str(),     CTag_Jet_CvsL     );
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsLN")) tree->SetBranchAddress((name+"CTag_Jet_CvsLN").c_str(),    CTag_Jet_CvsLN    );
      if(variableParser.isToBeStored(name+"CTag_Jet_CvsLP")) tree->SetBranchAddress((name+"CTag_Jet_CvsLP").c_str(),    CTag_Jet_CvsLP    );

      if(variableParser.isToBeStored(name+"Jet_nFirstTrkCTagVar")) tree->SetBranchAddress((name+"Jet_nFirstTrkCTagVar").c_str()        ,Jet_nFirstTrkCTagVar               );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkCTagVar")) tree->SetBranchAddress((name+"Jet_nLastTrkCTagVar").c_str()         ,Jet_nLastTrkCTagVar                     );
      if(variableParser.isToBeStored(name+"Jet_nFirstTrkEtaRelCTagVar")) tree->SetBranchAddress((name+"Jet_nFirstTrkEtaRelCTagVar").c_str()        ,Jet_nFirstTrkEtaRelCTagVar               );
      if(variableParser.isToBeStored(name+"Jet_nLastTrkEtaRelCTagVar")) tree->SetBranchAddress((name+"Jet_nLastTrkEtaRelCTagVar").c_str()        ,Jet_nLastTrkEtaRelCTagVar               );
      if(variableParser.isToBeStored(name+"Jet_nFirstLepCTagVar")) tree->SetBranchAddress((name+"Jet_nFirstLepCTagVar").c_str()        ,Jet_nFirstLepCTagVar               );
      if(variableParser.isToBeStored(name+"Jet_nLastLepCTagVar")) tree->SetBranchAddress((name+"Jet_nLastLepCTagVar").c_str()         ,Jet_nLastLepCTagVar                     );

      if(variableParser.isToBeStored(name+"CTag_jetNTracks")) tree->SetBranchAddress((name+"CTag_jetNTracks").c_str()               ,CTag_jetNTracks                             );
      if(variableParser.isToBeStored(name+"CTag_jetNTracksEtaRel")) tree->SetBranchAddress((name+"CTag_jetNTracksEtaRel").c_str()         ,CTag_jetNTracksEtaRel                 );
      if(variableParser.isToBeStored(name+"CTag_trackSumJetEtRatio")) tree->SetBranchAddress((name+"CTag_trackSumJetEtRatio").c_str()       ,CTag_trackSumJetEtRatio             );
      if(variableParser.isToBeStored(name+"CTag_trackSumJetDeltaR")) tree->SetBranchAddress((name+"CTag_trackSumJetDeltaR").c_str()        ,CTag_trackSumJetDeltaR               );
      if(variableParser.isToBeStored(name+"CTag_trackSip2dSigAboveCharm")) tree->SetBranchAddress((name+"CTag_trackSip2dSigAboveCharm").c_str()  ,CTag_trackSip2dSigAboveCharm   );
      if(variableParser.isToBeStored(name+"CTag_trackSip3dSigAboveCharm")) tree->SetBranchAddress((name+"CTag_trackSip3dSigAboveCharm").c_str()  ,CTag_trackSip3dSigAboveCharm   );
      if(variableParser.isToBeStored(name+"CTag_vertexCategory")) tree->SetBranchAddress((name+"CTag_vertexCategory").c_str()           ,CTag_vertexCategory                     );
      if(variableParser.isToBeStored(name+"CTag_jetNSecondaryVertices")) tree->SetBranchAddress((name+"CTag_jetNSecondaryVertices").c_str()    ,CTag_jetNSecondaryVertices       );
      if(variableParser.isToBeStored(name+"CTag_vertexMass")) tree->SetBranchAddress((name+"CTag_vertexMass").c_str()               ,CTag_vertexMass                             );
      if(variableParser.isToBeStored(name+"CTag_vertexNTracks")) tree->SetBranchAddress((name+"CTag_vertexNTracks").c_str()            ,CTag_vertexNTracks                       );
      if(variableParser.isToBeStored(name+"CTag_vertexEnergyRatio")) tree->SetBranchAddress((name+"CTag_vertexEnergyRatio").c_str()        ,CTag_vertexEnergyRatio               );
      if(variableParser.isToBeStored(name+"CTag_vertexJetDeltaR")) tree->SetBranchAddress((name+"CTag_vertexJetDeltaR").c_str()          ,CTag_vertexJetDeltaR                   );
      if(variableParser.isToBeStored(name+"CTag_flightDistance2dSig")) tree->SetBranchAddress((name+"CTag_flightDistance2dSig").c_str()      ,CTag_flightDistance2dSig           );
      if(variableParser.isToBeStored(name+"CTag_flightDistance3dSig")) tree->SetBranchAddress((name+"CTag_flightDistance3dSig").c_str()      ,CTag_flightDistance3dSig           );
      if(variableParser.isToBeStored(name+"CTag_massVertexEnergyFraction")) tree->SetBranchAddress((name+"CTag_massVertexEnergyFraction").c_str(),     CTag_massVertexEnergyFraction     );
      if(variableParser.isToBeStored(name+"CTag_vertexBoostOverSqrtJetPt")) tree->SetBranchAddress((name+"CTag_vertexBoostOverSqrtJetPt").c_str(),     CTag_vertexBoostOverSqrtJetPt     );

      if(variableParser.isToBeStored(name+"nTrkCTagVar")) tree->SetBranchAddress((name+"nTrkCTagVar").c_str()               ,&nTrkCTagVar                                                    );
      if(variableParser.isToBeStored(name+"nTrkEtaRelCTagVar")) tree->SetBranchAddress((name+"nTrkEtaRelCTagVar").c_str()         ,&nTrkEtaRelCTagVar                                        );
      if(variableParser.isToBeStored(name+"CTag_trackPtRel")) tree->SetBranchAddress((name+"CTag_trackPtRel").c_str()        ,CTag_trackPtRel                );
      if(variableParser.isToBeStored(name+"CTag_trackPPar")) tree->SetBranchAddress((name+"CTag_trackPPar").c_str()         ,CTag_trackPPar                  );
      if(variableParser.isToBeStored(name+"CTag_trackDeltaR")) tree->SetBranchAddress((name+"CTag_trackDeltaR").c_str()       ,CTag_trackDeltaR              );
      if(variableParser.isToBeStored(name+"CTag_trackPtRatio")) tree->SetBranchAddress((name+"CTag_trackPtRatio").c_str()      ,CTag_trackPtRatio            );
      if(variableParser.isToBeStored(name+"CTag_trackPParRatio")) tree->SetBranchAddress((name+"CTag_trackPParRatio").c_str()    ,CTag_trackPParRatio        );
      if(variableParser.isToBeStored(name+"CTag_trackSip2dSig")) tree->SetBranchAddress((name+"CTag_trackSip2dSig").c_str()     ,CTag_trackSip2dSig          );
      if(variableParser.isToBeStored(name+"CTag_trackSip3dSig")) tree->SetBranchAddress((name+"CTag_trackSip3dSig").c_str()     ,CTag_trackSip3dSig          );
      if(variableParser.isToBeStored(name+"CTag_trackDecayLenVal")) tree->SetBranchAddress((name+"CTag_trackDecayLenVal").c_str()  ,CTag_trackDecayLenVal    );
      if(variableParser.isToBeStored(name+"CTag_trackJetDistVal")) tree->SetBranchAddress((name+"CTag_trackJetDistVal").c_str()   ,CTag_trackJetDistVal      );
      if(variableParser.isToBeStored(name+"CTag_trackEtaRel")) tree->SetBranchAddress((name+"CTag_trackEtaRel").c_str()       ,CTag_trackEtaRel        );

      if(variableParser.isToBeStored(name+"CTag_vertexLeptonCategory")) tree->SetBranchAddress((name+"CTag_vertexLeptonCategory").c_str()    ,CTag_vertexLeptonCategory );
      if(variableParser.isToBeStored(name+"CTag_jetNLeptons")) tree->SetBranchAddress((name+"CTag_jetNLeptons").c_str()        ,CTag_jetNLeptons               );
      if(variableParser.isToBeStored(name+"nLeptons")) tree->SetBranchAddress((name+"nLeptons").c_str()               ,&nLeptons              );
      if(variableParser.isToBeStored(name+"CTag_leptonPtRel")) tree->SetBranchAddress((name+"CTag_leptonPtRel").c_str()         ,CTag_leptonPtRel         );
      if(variableParser.isToBeStored(name+"CTag_leptonSip3d")) tree->SetBranchAddress((name+"CTag_leptonSip3d").c_str()         ,CTag_leptonSip3d         );
      if(variableParser.isToBeStored(name+"CTag_leptonDeltaR")) tree->SetBranchAddress((name+"CTag_leptonDeltaR").c_str()         ,CTag_leptonDeltaR         );
      if(variableParser.isToBeStored(name+"CTag_leptonRatioRel")) tree->SetBranchAddress((name+"CTag_leptonRatioRel").c_str()         ,CTag_leptonRatioRel         );
      if(variableParser.isToBeStored(name+"CTag_leptonEtaRel")) tree->SetBranchAddress((name+"CTag_leptonEtaRel").c_str()         ,CTag_leptonEtaRel         );
      if(variableParser.isToBeStored(name+"CTag_leptonRatio")) tree->SetBranchAddress((name+"CTag_leptonRatio").c_str()         ,CTag_leptonRatio         );

      if(variableParser.isToBeStored(name+"Jet_FatJetIdx")) tree->SetBranchAddress((name+"Jet_FatJetIdx").c_str(),   Jet_FatJetIdx   );

      if(variableParser.isToBeStored(name+"Jet_ptSoftDrop")) tree->SetBranchAddress((name+"Jet_ptSoftDrop").c_str(),    Jet_ptSoftDrop    );
      if(variableParser.isToBeStored(name+"Jet_etaSoftDrop")) tree->SetBranchAddress((name+"Jet_etaSoftDrop").c_str(),   Jet_etaSoftDrop   );
      if(variableParser.isToBeStored(name+"Jet_phiSoftDrop")) tree->SetBranchAddress((name+"Jet_phiSoftDrop").c_str(),   Jet_phiSoftDrop   );
      if(variableParser.isToBeStored(name+"Jet_massSoftDrop")) tree->SetBranchAddress((name+"Jet_massSoftDrop").c_str(),  Jet_massSoftDrop  );
      if(variableParser.isToBeStored(name+"Jet_jecF0SoftDrop")) tree->SetBranchAddress((name+"Jet_jecF0SoftDrop").c_str(), Jet_jecF0SoftDrop );
      if(variableParser.isToBeStored(name+"Jet_ptPruned")) tree->SetBranchAddress((name+"Jet_ptPruned").c_str(),    Jet_ptPruned    );
      if(variableParser.isToBeStored(name+"Jet_etaPruned")) tree->SetBranchAddress((name+"Jet_etaPruned").c_str(),   Jet_etaPruned   );
      if(variableParser.isToBeStored(name+"Jet_phiPruned")) tree->SetBranchAddress((name+"Jet_phiPruned").c_str(),   Jet_phiPruned   );
      if(variableParser.isToBeStored(name+"Jet_massPruned")) tree->SetBranchAddress((name+"Jet_massPruned").c_str(),  Jet_massPruned  );
      if(variableParser.isToBeStored(name+"Jet_jecF0Pruned")) tree->SetBranchAddress((name+"Jet_jecF0Pruned").c_str(), Jet_jecF0Pruned );
      if(variableParser.isToBeStored(name+"Jet_tau1")) tree->SetBranchAddress((name+"Jet_tau1").c_str(),        Jet_tau1        );
      if(variableParser.isToBeStored(name+"Jet_tau2")) tree->SetBranchAddress((name+"Jet_tau2").c_str(),        Jet_tau2        );
      if(variableParser.isToBeStored(name+"Jet_tauAxis1_px")) tree->SetBranchAddress((name+"Jet_tauAxis1_px").c_str(), Jet_tauAxis1_px );
      if(variableParser.isToBeStored(name+"Jet_tauAxis1_py")) tree->SetBranchAddress((name+"Jet_tauAxis1_py").c_str(), Jet_tauAxis1_py );
      if(variableParser.isToBeStored(name+"Jet_tauAxis1_pz")) tree->SetBranchAddress((name+"Jet_tauAxis1_pz").c_str(), Jet_tauAxis1_pz );
      if(variableParser.isToBeStored(name+"Jet_tauAxis2_px")) tree->SetBranchAddress((name+"Jet_tauAxis2_px").c_str(), Jet_tauAxis2_px );
      if(variableParser.isToBeStored(name+"Jet_tauAxis2_py")) tree->SetBranchAddress((name+"Jet_tauAxis2_py").c_str(), Jet_tauAxis2_py );
      if(variableParser.isToBeStored(name+"Jet_tauAxis2_pz")) tree->SetBranchAddress((name+"Jet_tauAxis2_pz").c_str(), Jet_tauAxis2_pz );
      if(variableParser.isToBeStored(name+"Jet_z_ratio")) tree->SetBranchAddress((name+"Jet_z_ratio").c_str(),          Jet_z_ratio                   );
      if(variableParser.isToBeStored(name+"Jet_nTracks_fat")) tree->SetBranchAddress((name+"Jet_nTracks_fat").c_str(),      Jet_nTracks_fat           );
      if(variableParser.isToBeStored(name+"Jet_nSV_fat")) tree->SetBranchAddress((name+"Jet_nSV_fat").c_str(),          Jet_nSV_fat                   );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackEtaRel_0")) tree->SetBranchAddress((name+"Jet_tau1_trackEtaRel_0").c_str(),           Jet_tau1_trackEtaRel_0                      );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackEtaRel_1")) tree->SetBranchAddress((name+"Jet_tau1_trackEtaRel_1").c_str(),           Jet_tau1_trackEtaRel_1                      );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackEtaRel_2")) tree->SetBranchAddress((name+"Jet_tau1_trackEtaRel_2").c_str(),           Jet_tau1_trackEtaRel_2                      );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackEtaRel_0")) tree->SetBranchAddress((name+"Jet_tau2_trackEtaRel_0").c_str(),           Jet_tau2_trackEtaRel_0                      );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackEtaRel_1")) tree->SetBranchAddress((name+"Jet_tau2_trackEtaRel_1").c_str(),           Jet_tau2_trackEtaRel_1                      );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackEtaRel_2")) tree->SetBranchAddress((name+"Jet_tau2_trackEtaRel_2").c_str(),           Jet_tau2_trackEtaRel_2                      );
      if(variableParser.isToBeStored(name+"Jet_tau1_nSecondaryVertices")) tree->SetBranchAddress((name+"Jet_tau1_nSecondaryVertices").c_str(),      Jet_tau1_nSecondaryVertices            );
      if(variableParser.isToBeStored(name+"Jet_tau2_nSecondaryVertices")) tree->SetBranchAddress((name+"Jet_tau2_nSecondaryVertices").c_str(),      Jet_tau2_nSecondaryVertices            );
      if(variableParser.isToBeStored(name+"Jet_tau1_flightDistance2dSig")) tree->SetBranchAddress((name+"Jet_tau1_flightDistance2dSig").c_str(),     Jet_tau1_flightDistance2dSig          );
      if(variableParser.isToBeStored(name+"Jet_tau2_flightDistance2dSig")) tree->SetBranchAddress((name+"Jet_tau2_flightDistance2dSig").c_str(),     Jet_tau2_flightDistance2dSig          );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexDeltaR")) tree->SetBranchAddress((name+"Jet_tau1_vertexDeltaR").c_str(),            Jet_tau1_vertexDeltaR                        );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexDeltaR")) tree->SetBranchAddress((name+"Jet_tau2_vertexDeltaR").c_str(),            Jet_tau2_vertexDeltaR                        );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexEnergyRatio")) tree->SetBranchAddress((name+"Jet_tau1_vertexEnergyRatio").c_str(),       Jet_tau1_vertexEnergyRatio              );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexEnergyRatio")) tree->SetBranchAddress((name+"Jet_tau2_vertexEnergyRatio").c_str(),       Jet_tau2_vertexEnergyRatio              );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexMass")) tree->SetBranchAddress((name+"Jet_tau1_vertexMass").c_str(),              Jet_tau1_vertexMass                            );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexMass")) tree->SetBranchAddress((name+"Jet_tau2_vertexMass").c_str(),              Jet_tau2_vertexMass                            );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexMass_corrected")) tree->SetBranchAddress((name+"Jet_tau1_vertexMass_corrected").c_str(),    Jet_tau1_vertexMass_corrected        );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexMass_corrected")) tree->SetBranchAddress((name+"Jet_tau2_vertexMass_corrected").c_str(),    Jet_tau2_vertexMass_corrected        );
      if(variableParser.isToBeStored(name+"Jet_tau1_vertexNTracks")) tree->SetBranchAddress((name+"Jet_tau1_vertexNTracks").c_str(),           Jet_tau1_vertexNTracks                      );
      if(variableParser.isToBeStored(name+"Jet_tau2_vertexNTracks")) tree->SetBranchAddress((name+"Jet_tau2_vertexNTracks").c_str(),           Jet_tau2_vertexNTracks                      );
      if(variableParser.isToBeStored(name+"Jet_DoubleSV")) tree->SetBranchAddress((name+"Jet_DoubleSV").c_str(),         Jet_DoubleSV         );
      if(variableParser.isToBeStored(name+"Jet_DeepDoubleB")) tree->SetBranchAddress((name+"Jet_DeepDoubleB").c_str(),      Jet_DeepDoubleB      );


      if(variableParser.isToBeStored(name+"Track_lengthTau")) tree->SetBranchAddress((name+"Track_lengthTau").c_str()     ,Track_lengthTau       );
      if(variableParser.isToBeStored(name+"Track_distTau")) tree->SetBranchAddress((name+"Track_distTau").c_str()       ,Track_distTau         );

      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_3")) tree->SetBranchAddress((name+"Jet_trackSip3dSig_3").c_str(),                Jet_trackSip3dSig_3                         );
      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_2")) tree->SetBranchAddress((name+"Jet_trackSip3dSig_2").c_str(),                Jet_trackSip3dSig_2                         );
      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_1")) tree->SetBranchAddress((name+"Jet_trackSip3dSig_1").c_str(),                Jet_trackSip3dSig_1                         );
      if(variableParser.isToBeStored(name+"Jet_trackSip3dSig_0")) tree->SetBranchAddress((name+"Jet_trackSip3dSig_0").c_str(),                Jet_trackSip3dSig_0                         );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveCharm_0")) tree->SetBranchAddress((name+"Jet_trackSip2dSigAboveCharm_0").c_str(),      Jet_trackSip2dSigAboveCharm_0                   );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveCharm_1")) tree->SetBranchAddress((name+"Jet_trackSip2dSigAboveCharm_1").c_str(),      Jet_trackSip2dSigAboveCharm_1                   );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveBottom_0")) tree->SetBranchAddress((name+"Jet_trackSip2dSigAboveBottom_0").c_str(),     Jet_trackSip2dSigAboveBottom_0                 );
      if(variableParser.isToBeStored(name+"Jet_trackSip2dSigAboveBottom_1")) tree->SetBranchAddress((name+"Jet_trackSip2dSigAboveBottom_1").c_str(),     Jet_trackSip2dSigAboveBottom_1                 );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackSip3dSig_0")) tree->SetBranchAddress((name+"Jet_tau1_trackSip3dSig_0").c_str(),           Jet_tau1_trackSip3dSig_0                    );
      if(variableParser.isToBeStored(name+"Jet_tau1_trackSip3dSig_1")) tree->SetBranchAddress((name+"Jet_tau1_trackSip3dSig_1").c_str(),           Jet_tau1_trackSip3dSig_1                    );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackSip3dSig_0")) tree->SetBranchAddress((name+"Jet_tau2_trackSip3dSig_0").c_str(),           Jet_tau2_trackSip3dSig_0                    );
      if(variableParser.isToBeStored(name+"Jet_tau2_trackSip3dSig_1")) tree->SetBranchAddress((name+"Jet_tau2_trackSip3dSig_1").c_str(),           Jet_tau2_trackSip3dSig_1                    );
    }

};


class SubJetInfoBranches {

  public :

    int   Jet_nSubJets[nMaxJets_];
    int   Jet_nFirstSJ[nMaxJets_];
    int   Jet_nLastSJ[nMaxJets_];
    int   Jet_nsharedtracks[nMaxJets_];
    int   Jet_nsubjettracks[nMaxJets_];
    int   Jet_nsharedsubjettracks[nMaxJets_];

    int   nSubJet;
    int   SubJetIdx[nMaxJets_];


    void RegisterBranches(TTree *tree, VariableParser variableParser, std::string name="", std::string postfix="") {
      if(name!="") name += ".";
      if(postfix!="") postfix = "_" + postfix;

      if(variableParser.isToBeStored(name+"Jet_nSubJets"+postfix)) tree->Branch((name+"Jet_nSubJets"+postfix).c_str(),            Jet_nSubJets            ,(name+"Jet_nSubJets"+postfix+"["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nFirstSJ"+postfix)) tree->Branch((name+"Jet_nFirstSJ"+postfix).c_str(),            Jet_nFirstSJ            ,(name+"Jet_nFirstSJ"+postfix+"["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nLastSJ"+postfix)) tree->Branch((name+"Jet_nLastSJ"+postfix).c_str(),             Jet_nLastSJ             ,(name+"Jet_nLastSJ"+postfix+"["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nsharedtracks"+postfix)) tree->Branch((name+"Jet_nsharedtracks"+postfix).c_str(),       Jet_nsharedtracks       ,(name+"Jet_nsharedtracks"+postfix+"["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nsubjettracks"+postfix)) tree->Branch((name+"Jet_nsubjettracks"+postfix).c_str(),       Jet_nsubjettracks       ,(name+"Jet_nsubjettracks"+postfix+"["+name+"nJet]/I").c_str());
      if(variableParser.isToBeStored(name+"Jet_nsharedsubjettracks"+postfix)) tree->Branch((name+"Jet_nsharedsubjettracks"+postfix).c_str(), Jet_nsharedsubjettracks ,(name+"Jet_nsharedsubjettracks"+postfix+"["+name+"nJet]/I").c_str());

      if(variableParser.isToBeStored(name+"nSubJet"+postfix)) tree->Branch((name+"nSubJet"+postfix).c_str(),                 &nSubJet                ,(name+"nSubJet"+postfix+"/I").c_str());
      if(variableParser.isToBeStored(name+"SubJetIdx"+postfix)) tree->Branch((name+"SubJetIdx"+postfix).c_str(),               SubJetIdx               ,(name+"SubJetIdx"+postfix+"["+name+"nSubJet"+postfix+"]/I").c_str());
    }

    void ReadBranches(TTree *tree, VariableParser variableParser, std::string name="", std::string postfix="") {
      if(name!="") name += ".";
      if(postfix!="") postfix = "_" + postfix;

      if(variableParser.isToBeStored(name+"Jet_nSubJets"+postfix)) tree->SetBranchAddress((name+"Jet_nSubJets"+postfix).c_str(),            Jet_nSubJets           );
      if(variableParser.isToBeStored(name+"Jet_nFirstSJ"+postfix)) tree->SetBranchAddress((name+"Jet_nFirstSJ"+postfix).c_str(),            Jet_nFirstSJ           );
      if(variableParser.isToBeStored(name+"Jet_nLastSJ"+postfix)) tree->SetBranchAddress((name+"Jet_nLastSJ"+postfix).c_str(),             Jet_nLastSJ             );
      if(variableParser.isToBeStored(name+"Jet_nsharedtracks"+postfix)) tree->SetBranchAddress((name+"Jet_nsharedtracks"+postfix).c_str(),       Jet_nsharedtracks );
      if(variableParser.isToBeStored(name+"Jet_nsubjettracks"+postfix)) tree->SetBranchAddress((name+"Jet_nsubjettracks"+postfix).c_str(),       Jet_nsubjettracks );
      if(variableParser.isToBeStored(name+"Jet_nsharedsubjettracks"+postfix)) tree->SetBranchAddress((name+"Jet_nsharedsubjettracks"+postfix).c_str(), Jet_nsharedsubjettracks );

      if(variableParser.isToBeStored(name+"nSubJet"+postfix)) tree->SetBranchAddress((name+"nSubJet"+postfix).c_str(),                 &nSubJet                );
      if(variableParser.isToBeStored(name+"SubJetIdx"+postfix)) tree->SetBranchAddress((name+"SubJetIdx"+postfix).c_str(),               SubJetIdx             );
    }

};

#endif
