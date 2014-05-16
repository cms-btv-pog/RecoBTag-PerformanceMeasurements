#ifndef JETINFOBRANCHES_H
#define JETINFOBRANCHES_H

#include <TTree.h>

const UInt_t nMaxJets_ = 10000;
const UInt_t nMaxTrk_  = 100000;
const UInt_t nMaxMuons_= 10000;
const UInt_t nMaxElectrons_= 10000;
const UInt_t nMaxSVs_= 10000;

class JetInfoBranches {

  public :

    int   nJet;
    float Jet_pt[nMaxJets_];
    float Jet_genpt[nMaxJets_];
    float Jet_residual[nMaxJets_];
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
    float Jet_ProbaN[nMaxJets_];
    float Jet_ProbaP[nMaxJets_];
    float Jet_Proba[nMaxJets_];
    float Jet_BprobN[nMaxJets_];
    float Jet_BprobP[nMaxJets_];
    float Jet_Bprob[nMaxJets_];
    float Jet_SvxN[nMaxJets_];
    float Jet_Svx[nMaxJets_];
    int   Jet_SvxNTracks[nMaxJets_];
    int   Jet_SvxTracks[nMaxJets_];
    float Jet_SvxNHP[nMaxJets_];
    float Jet_SvxHP[nMaxJets_];
    float Jet_SvxMass[nMaxJets_];
    float Jet_CombSvxN[nMaxJets_];
    float Jet_CombSvxP[nMaxJets_];
    float Jet_CombSvx[nMaxJets_];
    float Jet_RetCombSvxN[nMaxJets_];
    float Jet_RetCombSvxP[nMaxJets_];
    float Jet_RetCombSvx[nMaxJets_];
    float Jet_CombCSVJP_N[nMaxJets_];
    float Jet_CombCSVJP_P[nMaxJets_];
    float Jet_CombCSVJP[nMaxJets_];
    float Jet_CombCSVSL_N[nMaxJets_];
    float Jet_CombCSVSL_P[nMaxJets_];
    float Jet_CombCSVSL[nMaxJets_];
    float Jet_CombCSVJPSL_N[nMaxJets_];
    float Jet_CombCSVJPSL_P[nMaxJets_];
    float Jet_CombCSVJPSL[nMaxJets_];
    float Jet_SimpIVF_HP[nMaxJets_];
    float Jet_SimpIVF_HE[nMaxJets_];
    float Jet_DoubIVF_HE[nMaxJets_];
    float Jet_CombIVF[nMaxJets_];
    float Jet_CombIVF_P[nMaxJets_];
    float Jet_SoftMuN[nMaxJets_];
    float Jet_SoftMuP[nMaxJets_];
    float Jet_SoftMu[nMaxJets_];
    float Jet_SoftElN[nMaxJets_];
    float Jet_SoftElP[nMaxJets_];
    float Jet_SoftEl[nMaxJets_];
    int   Jet_hist1[nMaxJets_];
    int   Jet_hist2[nMaxJets_];
    int   Jet_hist3[nMaxJets_];
    int   Jet_histJet[nMaxJets_];
    int   Jet_histSvx[nMaxJets_];
    int   Jet_ntracks[nMaxJets_];
    int   Jet_nseltracks[nMaxJets_];
    int   Jet_nsubjettracks[nMaxJets_];
    int   Jet_nsharedsubjettracks[nMaxJets_];
    int   Jet_nsharedtracks[nMaxJets_];
    int   Jet_flavour[nMaxJets_];
    int   Jet_nFirstTrack[nMaxJets_];
    int   Jet_nLastTrack[nMaxJets_];
    int   Jet_nFirstSV[nMaxJets_];
    int   Jet_nLastSV[nMaxJets_];
    int   Jet_nFirstTrkInc[nMaxJets_];
    int   Jet_nLastTrkInc[nMaxJets_];
    int   Jet_SV_multi[nMaxJets_];
    int   Jet_VtxCat[nMaxJets_];
    int   Jet_looseID[nMaxJets_];
    int   Jet_tightID[nMaxJets_];
    int   Jet_FatJetIdx[nMaxJets_];
    float Jet_ptGroomed[nMaxJets_];
    float Jet_jesGroomed[nMaxJets_];
    float Jet_etaGroomed[nMaxJets_];
    float Jet_phiGroomed[nMaxJets_];
    float Jet_massGroomed[nMaxJets_];
    float Jet_tau1[nMaxJets_];
    float Jet_tau2[nMaxJets_];
    int   Jet_nSubJets[nMaxJets_];
    int   Jet_nFirstSJ[nMaxJets_];
    int   Jet_nLastSJ[nMaxJets_];
    int   Jet_nFirstTrkTagVar[nMaxJets_];
    int   Jet_nLastTrkTagVar[nMaxJets_];
    int   Jet_nFirstSVTagVar[nMaxJets_];
    int   Jet_nLastSVTagVar[nMaxJets_];
    int   Jet_nFirstTrkTagVarCSV[nMaxJets_];
    int   Jet_nLastTrkTagVarCSV[nMaxJets_];
    int   Jet_nFirstTrkEtaRelTagVarCSV[nMaxJets_];
    int   Jet_nLastTrkEtaRelTagVarCSV[nMaxJets_];

    int   nSubJet;
    int   SubJetIdx[nMaxJets_];

    int   nTrack;
    float Track_dxy[nMaxTrk_];
    float Track_dz[nMaxTrk_];
    float Track_zIP[nMaxTrk_];
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
    int   Track_category[nMaxTrk_];

    int   nTrkInc;
    float TrkInc_pt[nMaxTrk_];
    float TrkInc_eta[nMaxTrk_];
    float TrkInc_phi[nMaxTrk_];
    float TrkInc_ptrel[nMaxTrk_];
    float TrkInc_IPsig[nMaxTrk_];
    float TrkInc_IP[nMaxTrk_];

    int   nMuon;
    int   Muon_IdxJet[nMaxMuons_];
    int   Muon_nMuHit[nMaxMuons_];
    int   Muon_nTkHit[nMaxMuons_];
    int   Muon_nPixHit[nMaxMuons_];
    int   Muon_nOutHit[nMaxMuons_];
    int   Muon_isGlobal[nMaxMuons_];
    int   Muon_nMatched[nMaxMuons_];
    float Muon_chi2[nMaxMuons_];
    float Muon_chi2Tk[nMaxMuons_];
    float Muon_pt[nMaxMuons_];
    float Muon_eta[nMaxMuons_];
    float Muon_phi[nMaxMuons_];
    float Muon_ptrel[nMaxMuons_];
    float Muon_vz[nMaxMuons_];
    int   Muon_hist[nMaxMuons_];
    int   Muon_TrackIdx[nMaxMuons_];
    float Muon_IPsig[nMaxMuons_];
    float Muon_IP[nMaxMuons_];
    float Muon_IP2Dsig[nMaxMuons_];
    float Muon_IP2D[nMaxMuons_];
    float Muon_Proba[nMaxMuons_];
    float Muon_deltaR[nMaxMuons_];
    float Muon_ratio[nMaxMuons_];
    float Muon_ratioRel[nMaxMuons_];

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
    int   PFMuon_IdxJet[nMaxElectrons_];
    float PFMuon_pt[nMaxElectrons_];
    float PFMuon_eta[nMaxElectrons_];
    float PFMuon_phi[nMaxElectrons_];
    float PFMuon_ptrel[nMaxElectrons_];
    float PFMuon_ratio[nMaxElectrons_];
    float PFMuon_ratioRel[nMaxElectrons_];
    float PFMuon_deltaR[nMaxElectrons_];
    float PFMuon_IP[nMaxElectrons_];
    float PFMuon_IP2D[nMaxElectrons_];
    int   PFMuon_GoodQuality[nMaxElectrons_];

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
    float SV_energy_ratio[nMaxSVs_];
    float SV_aboveC[nMaxSVs_];
    float SV_vtx_pt[nMaxSVs_];
    float SV_flight2D[nMaxSVs_];
    float SV_flight2DErr[nMaxSVs_];
    float SV_totCharge[nMaxSVs_];
    float SV_vtxDistJetAxis[nMaxSVs_];
    int   SV_nTrk[nMaxSVs_];
    int   SV_nTrk_firstVxt[nMaxSVs_];
    float SV_mass[nMaxSVs_];
    float SV_vtx_eta[nMaxSVs_];
    float SV_vtx_phi[nMaxSVs_];

    // TagInfo TaggingVariables
    // per jet
    float TagVar_jetNTracks[nMaxJets_];                              // tracks associated to jet
    float TagVar_jetNSecondaryVertices[nMaxJets_];                   // number of reconstructed possible secondary vertices in jet
    float TagVar_chargedHadronEnergyFraction[nMaxJets_];             // fraction of the jet energy coming from charged hadrons
    float TagVar_neutralHadronEnergyFraction[nMaxJets_];             // fraction of the jet energy coming from neutral hadrons
    float TagVar_photonEnergyFraction[nMaxJets_];                    // fraction of the jet energy coming from photons
    float TagVar_electronEnergyFraction[nMaxJets_];                  // fraction of the jet energy coming from electrons
    float TagVar_muonEnergyFraction[nMaxJets_];                      // fraction of the jet energy coming from muons
    float TagVar_chargedHadronMultiplicity[nMaxJets_];               // number of charged hadrons in the jet
    float TagVar_neutralHadronMultiplicity[nMaxJets_];               // number of neutral hadrons in the jet
    float TagVar_photonMultiplicity[nMaxJets_];                      // number of photons in the jet
    float TagVar_electronMultiplicity[nMaxJets_];                    // number of electrons in the jet
    float TagVar_muonMultiplicity[nMaxJets_];                        // number of muons in the jet
    // per jet per track
    int   nTrkTagVar;
    float TagVar_trackMomentum[nMaxTrk_];                            // track momentum
    float TagVar_trackEta[nMaxTrk_];                                 // track pseudorapidity
    float TagVar_trackPhi[nMaxTrk_];                                 // track polar angle
    float TagVar_trackPtRel[nMaxTrk_];                               // track transverse momentum, relative to the jet axis
    float TagVar_trackPPar[nMaxTrk_];                                // track parallel momentum, along the jet axis
    float TagVar_trackEtaRel[nMaxTrk_];                              // track pseudorapidity, relative to the jet axis
    float TagVar_trackDeltaR[nMaxTrk_];                              // track pseudoangular distance from the jet axis
    float TagVar_trackPtRatio[nMaxTrk_];                             // track transverse momentum, relative to the jet axis, normalized to its energy
    float TagVar_trackPParRatio[nMaxTrk_];                           // track parallel momentum, along the jet axis, normalized to its energy
    float TagVar_trackSip2dVal[nMaxTrk_];                            // track 2D signed impact parameter
    float TagVar_trackSip2dSig[nMaxTrk_];                            // track 2D signed impact parameter significance
    float TagVar_trackSip3dVal[nMaxTrk_];                            // track 3D signed impact parameter
    float TagVar_trackSip3dSig[nMaxTrk_];                            // track 3D signed impact parameter significance
    float TagVar_trackDecayLenVal[nMaxTrk_];                         // track decay length
    float TagVar_trackDecayLenSig[nMaxTrk_];                         // track decay length significance
    float TagVar_trackJetDistVal[nMaxTrk_];                          // minimum track approach distance to jet axis
    float TagVar_trackJetDistSig[nMaxTrk_];                          // minimum track approach distance to jet axis significance
    float TagVar_trackChi2[nMaxTrk_];                                // track fit chi2
    float TagVar_trackNTotalHits[nMaxTrk_];                          // number of valid total hits
    float TagVar_trackNPixelHits[nMaxTrk_];                          // number of valid pixel hits
    // per jet per secondary vertex
    int   nSVTagVar;
    float TagVar_vertexMass[nMaxSVs_];                               // mass of track sum at secondary vertex
    float TagVar_vertexNTracks[nMaxSVs_];                            // number of tracks at secondary vertex
    float TagVar_vertexJetDeltaR[nMaxSVs_];                          // pseudoangular distance between jet axis and secondary vertex direction
    float TagVar_flightDistance2dVal[nMaxSVs_];                      // transverse distance between primary and secondary vertex
    float TagVar_flightDistance2dSig[nMaxSVs_];                      // transverse distance significance between primary and secondary vertex
    float TagVar_flightDistance3dVal[nMaxSVs_];                      // distance between primary and secondary vertex
    float TagVar_flightDistance3dSig[nMaxSVs_];                      // distance significance between primary and secondary vertex

    // CSV TaggingVariables
    // per jet
    float TagVarCSV_trackJetPt[nMaxJets_];                           // track-based jet transverse momentum
    float TagVarCSV_jetNTracks[nMaxJets_];                           // tracks associated to jet
    float TagVarCSV_jetNTracksEtaRel[nMaxJets_];                     // tracks associated to jet for which trackEtaRel is calculated
    float TagVarCSV_trackSumJetEtRatio[nMaxJets_];                   // ratio of track sum transverse energy over jet energy
    float TagVarCSV_trackSumJetDeltaR[nMaxJets_];                    // pseudoangular distance between jet axis and track fourvector sum
    float TagVarCSV_trackSip2dValAboveCharm[nMaxJets_];              // track 2D signed impact parameter of first track lifting mass above charm
    float TagVarCSV_trackSip2dSigAboveCharm[nMaxJets_];              // track 2D signed impact parameter significance of first track lifting mass above charm
    float TagVarCSV_trackSip3dValAboveCharm[nMaxJets_];              // track 3D signed impact parameter of first track lifting mass above charm
    float TagVarCSV_trackSip3dSigAboveCharm[nMaxJets_];              // track 3D signed impact parameter significance of first track lifting mass above charm
    float TagVarCSV_vertexCategory[nMaxJets_];                       // category of secondary vertex (Reco, Pseudo, No)
    float TagVarCSV_jetNSecondaryVertices[nMaxJets_];                // number of reconstructed possible secondary vertices in jet
    float TagVarCSV_vertexMass[nMaxJets_];                           // mass of track sum at secondary vertex
    float TagVarCSV_vertexNTracks[nMaxJets_];                        // number of tracks at secondary vertex
    float TagVarCSV_vertexEnergyRatio[nMaxJets_];                    // ratio of energy at secondary vertex over total energy
    float TagVarCSV_vertexJetDeltaR[nMaxJets_];                      // pseudoangular distance between jet axis and secondary vertex direction
    float TagVarCSV_flightDistance2dVal[nMaxJets_];                  // transverse distance between primary and secondary vertex
    float TagVarCSV_flightDistance2dSig[nMaxJets_];                  // transverse distance significance between primary and secondary vertex
    float TagVarCSV_flightDistance3dVal[nMaxJets_];                  // distance between primary and secondary vertex
    float TagVarCSV_flightDistance3dSig[nMaxJets_];                  // distance significance between primary and secondary vertex
    // per jet per track
    int   nTrkTagVarCSV;
    int   nTrkEtaRelTagVarCSV;
    float TagVarCSV_trackMomentum[nMaxTrk_];                         // track momentum
    float TagVarCSV_trackEta[nMaxTrk_];                              // track pseudorapidity
    float TagVarCSV_trackPhi[nMaxTrk_];                              // track polar angle
    float TagVarCSV_trackPtRel[nMaxTrk_];                            // track transverse momentum, relative to the jet axis
    float TagVarCSV_trackPPar[nMaxTrk_];                             // track parallel momentum, along the jet axis
    float TagVarCSV_trackDeltaR[nMaxTrk_];                           // track pseudoangular distance from the jet axis
    float TagVarCSV_trackPtRatio[nMaxTrk_];                          // track transverse momentum, relative to the jet axis, normalized to its energy
    float TagVarCSV_trackPParRatio[nMaxTrk_];                        // track parallel momentum, along the jet axis, normalized to its energy
    float TagVarCSV_trackSip2dVal[nMaxTrk_];                         // track 2D signed impact parameter
    float TagVarCSV_trackSip2dSig[nMaxTrk_];                         // track 2D signed impact parameter significance
    float TagVarCSV_trackSip3dVal[nMaxTrk_];                         // track 3D signed impact parameter
    float TagVarCSV_trackSip3dSig[nMaxTrk_];                         // track 3D signed impact parameter significance
    float TagVarCSV_trackDecayLenVal[nMaxTrk_];                      // track decay length
    float TagVarCSV_trackDecayLenSig[nMaxTrk_];                      // track decay length significance
    float TagVarCSV_trackJetDistVal[nMaxTrk_];                       // minimum track approach distance to jet axis
    float TagVarCSV_trackJetDistSig[nMaxTrk_];                       // minimum track approach distance to jet axis significance
    float TagVarCSV_trackEtaRel[nMaxTrk_];                           // track pseudorapidity, relative to the jet axis


    void RegisterTree(TTree *tree, std::string name="") {
      if(name!="") name += ".";
      tree->Branch((name+"nJet").c_str(),            &nJet           ,(name+"nJet/I").c_str());
      tree->Branch((name+"Jet_pt").c_str(),          Jet_pt        ,(name+"Jet_pt["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_genpt").c_str(),       Jet_genpt       ,(name+"Jet_genpt["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_residual").c_str(),    Jet_residual    ,(name+"Jet_residual["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_jes").c_str(),         Jet_jes         ,(name+"Jet_jes["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_eta").c_str(),         Jet_eta         ,(name+"Jet_eta["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_phi").c_str(),         Jet_phi         ,(name+"Jet_phi["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_mass").c_str(),        Jet_mass        ,(name+"Jet_mass["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_ntracks").c_str(),     Jet_ntracks     ,(name+"Jet_ntracks["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nseltracks").c_str(),  Jet_nseltracks  ,(name+"Jet_nseltracks["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_flavour").c_str(),     Jet_flavour     ,(name+"Jet_flavour["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_Ip2N").c_str(),        Jet_Ip2N        ,(name+"Jet_Ip2N["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_Ip2P").c_str(),        Jet_Ip2P        ,(name+"Jet_Ip2P["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_Ip3N").c_str(),        Jet_Ip3N        ,(name+"Jet_Ip3N["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_Ip3P").c_str(),        Jet_Ip3P        ,(name+"Jet_Ip3P["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_ProbaN").c_str(),      Jet_ProbaN     ,(name+"Jet_ProbaN["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_ProbaP").c_str(),      Jet_ProbaP     ,(name+"Jet_ProbaP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_Proba").c_str(),       Jet_Proba      ,(name+"Jet_Proba["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_BprobN").c_str(),      Jet_BprobN     ,(name+"Jet_BprobN["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_BprobP").c_str(),      Jet_BprobP     ,(name+"Jet_BprobP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_Bprob").c_str(),       Jet_Bprob      ,(name+"Jet_Bprob["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SvxN").c_str(),        Jet_SvxN       ,(name+"Jet_SvxN["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_Svx").c_str(),         Jet_Svx        ,(name+"Jet_Svx["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SvxNHP").c_str(),      Jet_SvxNHP     ,(name+"Jet_SvxNHP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SvxHP").c_str(),       Jet_SvxHP      ,(name+"Jet_SvxHP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SvxMass").c_str(),     Jet_SvxMass    ,(name+"Jet_SvxMass["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombSvxN").c_str(),    Jet_CombSvxN   ,(name+"Jet_CombSvxN["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombSvxP").c_str(),    Jet_CombSvxP   ,(name+"Jet_CombSvxP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombSvx").c_str(),     Jet_CombSvx    ,(name+"Jet_CombSvx["+name+"nJet]/F").c_str());

      tree->Branch((name+"Jet_RetCombSvxN").c_str(), Jet_RetCombSvxN   ,(name+"Jet_RetCombSvxN["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_RetCombSvxP").c_str(), Jet_RetCombSvxP   ,(name+"Jet_RetCombSvxP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_RetCombSvx").c_str(), Jet_RetCombSvx    ,(name+"Jet_RetCombSvx["+name+"nJet]/F").c_str());

      tree->Branch((name+"Jet_CombCSVJP_N").c_str(), Jet_CombCSVJP_N   ,(name+"Jet_CombCSVJP_N["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombCSVJP_P").c_str(), Jet_CombCSVJP_P   ,(name+"Jet_CombCSVJP_P["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombCSVJP").c_str(),Jet_CombCSVJP      ,(name+"Jet_CombCSVJP["+name+"nJet]/F").c_str());

      tree->Branch((name+"Jet_CombCSVSL_N").c_str(), Jet_CombCSVSL_N    ,(name+"Jet_CombCSVSL_N["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombCSVSL_P").c_str(), Jet_CombCSVSL_P    ,(name+"Jet_CombCSVSL_P["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombCSVSL").c_str(),   Jet_CombCSVSL ,(name+"Jet_CombCSVSL["+name+"nJet]/F").c_str());

      tree->Branch((name+"Jet_CombCSVJPSL_N").c_str(), Jet_CombCSVJPSL_N  ,(name+"Jet_CombCSVJPSL_N["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombCSVJPSL_P").c_str(), Jet_CombCSVJPSL_P  ,(name+"Jet_CombCSVJPSL_P["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_CombCSVJPSL").c_str(), Jet_CombCSVJPSL    ,(name+"Jet_CombCSVJPSL["+name+"nJet]/F").c_str());

      // tree->Branch((name+"Jet_SimpIVF_HP").c_str(),  Jet_SimpIVF_HP  ,(name+"Jet_SimpIVF_HP["+name+"nJet]/F").c_str());
      // tree->Branch((name+"Jet_SimpIVF_HE").c_str(),  Jet_SimpIVF_HE  ,(name+"Jet_SimpIVF_HE["+name+"nJet]/F").c_str());
      // tree->Branch((name+"Jet_DoubIVF_HE").c_str(),  Jet_DoubIVF_HE  ,(name+"Jet_DoubIVF_HE["+name+"nJet]/F").c_str());
      // tree->Branch((name+"Jet_CombIVF").c_str(),     Jet_CombIVF     ,(name+"Jet_CombIVF["+name+"nJet]/F").c_str());
      // tree->Branch((name+"Jet_CombIVF_P").c_str(), Jet_CombIVF_P   ,(name+"Jet_CombIVF_P["+name+"nJet]/F").c_str());

      tree->Branch((name+"Jet_SoftMuN").c_str(),     Jet_SoftMuN     ,(name+"Jet_SoftMuN["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SoftMuP").c_str(),     Jet_SoftMuP     ,(name+"Jet_SoftMuP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SoftMu").c_str(),      Jet_SoftMu      ,(name+"Jet_SoftMu["+name+"nJet]/F").c_str());

      tree->Branch((name+"Jet_SoftElN").c_str(),     Jet_SoftElN     ,(name+"Jet_SoftElN["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SoftElP").c_str(),     Jet_SoftElP     ,(name+"Jet_SoftElP["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_SoftEl").c_str(),      Jet_SoftEl      ,(name+"Jet_SoftEl["+name+"nJet]/F").c_str());

      tree->Branch((name+"Jet_hist1").c_str(),       Jet_hist1       ,(name+"Jet_hist1["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_hist2").c_str(),       Jet_hist2       ,(name+"Jet_hist2["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_hist3").c_str(),       Jet_hist3       ,(name+"Jet_hist3["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_histJet").c_str(),     Jet_histJet     ,(name+"Jet_histJet["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_histSvx").c_str(),     Jet_histSvx     ,(name+"Jet_histSvx["+name+"nJet]/I").c_str());

      tree->Branch((name+"Jet_nFirstTrack").c_str(), Jet_nFirstTrack ,(name+"Jet_nFirstTrack["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nLastTrack").c_str(),  Jet_nLastTrack  ,(name+"Jet_nLastTrack["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nFirstSV").c_str(),    Jet_nFirstSV    ,(name+"Jet_nFirstSV["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nLastSV").c_str(),     Jet_nLastSV     ,(name+"Jet_nLastSV["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_SV_multi").c_str(),    Jet_SV_multi      ,(name+"Jet_SV_multi["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nFirstTrkInc").c_str(), Jet_nFirstTrkInc ,(name+"Jet_nFirstTrkInc["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nLastTrkInc").c_str(),  Jet_nLastTrkInc  ,(name+"Jet_nLastTrkInc["+name+"nJet]/I").c_str());

      tree->Branch((name+"Jet_VtxCat").c_str(),      Jet_VtxCat    ,(name+"Jet_VtxCat["+name+"nJet]/I").c_str() );
      tree->Branch((name+"Jet_looseID").c_str(),      Jet_looseID  ,(name+"Jet_looseID["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_tightID").c_str(),      Jet_tightID  ,(name+"Jet_tightID["+name+"nJet]/I").c_str());

      //--------------------------------------
      // Inclusive Track information for PtRel template
      //--------------------------------------
      tree->Branch((name+"nTrkInc").c_str()      ,&nTrkInc     ,(name+"nTrkInc/I").c_str());
      tree->Branch((name+"TrkInc_pt").c_str()    ,TrkInc_pt    ,(name+"TrkInc_pt["+name+"nTrkInc]/F").c_str());
      tree->Branch((name+"TrkInc_eta").c_str()   ,TrkInc_eta   ,(name+"TrkInc_eta["+name+"nTrkInc]/F").c_str());
      tree->Branch((name+"TrkInc_phi").c_str()   ,TrkInc_phi   ,(name+"TrkInc_phi["+name+"nTrkInc]/F").c_str());
      tree->Branch((name+"TrkInc_ptrel").c_str() ,TrkInc_ptrel ,(name+"TrkInc_ptrel["+name+"nTrkInc]/F").c_str());
      tree->Branch((name+"TrkInc_IPsig").c_str() ,TrkInc_IPsig ,(name+"TrkInc_IPsig["+name+"nTrkInc]/F").c_str());
      tree->Branch((name+"TrkInc_IP").c_str()    ,TrkInc_IP    ,(name+"TrkInc_IP["+name+"nTrkInc]/F").c_str());

      //--------------------------------------
      // muon information
      //--------------------------------------
      tree->Branch((name+"nMuon").c_str()         ,&nMuon        ,(name+"nMuon/I").c_str());
      tree->Branch((name+"Muon_IdxJet").c_str()   ,Muon_IdxJet   ,(name+"Muon_IdxJet["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_nMuHit").c_str()   ,Muon_nMuHit   ,(name+"Muon_nMuHit["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_nTkHit").c_str()   ,Muon_nTkHit   ,(name+"Muon_nTkHit["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_nPixHit").c_str()  ,Muon_nPixHit  ,(name+"Muon_nPixHit["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_nOutHit").c_str()  ,Muon_nOutHit  ,(name+"Muon_nOutHit["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_isGlobal").c_str() ,Muon_isGlobal ,(name+"Muon_isGlobal["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_nMatched").c_str() ,Muon_nMatched ,(name+"Muon_nMatched["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_chi2").c_str()	  ,Muon_chi2	 ,(name+"Muon_chi2["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_chi2Tk").c_str()   ,Muon_chi2Tk   ,(name+"Muon_chi2Tk["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_pt").c_str()       ,Muon_pt       ,(name+"Muon_pt["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_eta").c_str()      ,Muon_eta      ,(name+"Muon_eta["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_phi").c_str()      ,Muon_phi      ,(name+"Muon_phi["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_ptrel").c_str()    ,Muon_ptrel    ,(name+"Muon_ptrel["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_vz").c_str()	  ,Muon_vz	 ,(name+"Muon_vz["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_hist").c_str()	  ,Muon_hist	 ,(name+"Muon_hist["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_TrackIdx").c_str() ,Muon_TrackIdx ,(name+"Muon_TrackIdx["+name+"nMuon]/I").c_str());
      tree->Branch((name+"Muon_IPsig").c_str()    ,Muon_IPsig	 ,(name+"Muon_IPsig["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_IP").c_str()       ,Muon_IP       ,(name+"Muon_IP["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_IP2Dsig").c_str()  ,Muon_IP2Dsig  ,(name+"Muon_IP2Dsig["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_IP2D").c_str()     ,Muon_IP2D     ,(name+"Muon_IP2D["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_Proba").c_str()    ,Muon_Proba    ,(name+"Muon_Proba["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_deltaR").c_str()   ,Muon_deltaR   ,(name+"Muon_deltaR["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_ratio").c_str()    ,Muon_ratio    ,(name+"Muon_ratio["+name+"nMuon]/F").c_str());
      tree->Branch((name+"Muon_ratioRel").c_str() ,Muon_ratioRel ,(name+"Muon_ratioRel["+name+"nMuon]/F").c_str());

      //--------------------------------------
      // pf electron information
      //--------------------------------------
      tree->Branch((name+"nPFElectron").c_str()         ,&nPFElectron   ,(name+"nPFElectron/I").c_str());
      tree->Branch((name+"PFElectron_IdxJet").c_str()   ,PFElectron_IdxJet   ,(name+"PFElectron_IdxJet["+name+"nPFElectron]/I").c_str());
      tree->Branch((name+"PFElectron_pt").c_str()       ,PFElectron_pt       ,(name+"PFElectron_pt["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_eta").c_str()      ,PFElectron_eta      ,(name+"PFElectron_eta["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_phi").c_str()      ,PFElectron_phi      ,(name+"PFElectron_phi["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_ptrel").c_str()    ,PFElectron_ptrel    ,(name+"PFElectron_ptrel["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_deltaR").c_str()   ,PFElectron_deltaR   ,(name+"PFElectron_deltaR["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_ratio").c_str()    ,PFElectron_ratio    ,(name+"PFElectron_ratio["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_ratioRel").c_str() ,PFElectron_ratioRel ,(name+"PFElectron_ratioRel["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_IP").c_str()       ,PFElectron_IP       ,(name+"PFElectron_IP["+name+"nPFElectron]/F").c_str());
      tree->Branch((name+"PFElectron_IP2D").c_str()     ,PFElectron_IP2D     ,(name+"PFElectron_IP2D["+name+"nPFElectron]/F").c_str());

      //--------------------------------------
      // pf muon information
      //--------------------------------------
      tree->Branch((name+"nPFMuon").c_str()            ,&nPFMuon            ,(name+"nPFMuon/I").c_str());
      tree->Branch((name+"PFMuon_IdxJet").c_str()      ,PFMuon_IdxJet       ,(name+"PFMuon_IdxJet["+name+"nPFMuon]/I").c_str());
      tree->Branch((name+"PFMuon_pt").c_str()          ,PFMuon_pt           ,(name+"PFMuon_pt["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_eta").c_str()         ,PFMuon_eta          ,(name+"PFMuon_eta["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_phi").c_str()         ,PFMuon_phi          ,(name+"PFMuon_phi["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_ptrel").c_str()       ,PFMuon_ptrel        ,(name+"PFMuon_ptrel["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_deltaR").c_str()      ,PFMuon_deltaR       ,(name+"PFMuon_deltaR["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_ratio").c_str()       ,PFMuon_ratio        ,(name+"PFMuon_ratio["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_ratioRel").c_str()    ,PFMuon_ratioRel     ,(name+"PFMuon_ratioRel["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_IP").c_str()          ,PFMuon_IP           ,(name+"PFMuon_IP["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_IP2D").c_str()        ,PFMuon_IP2D         ,(name+"PFMuon_IP2D["+name+"nPFMuon]/F").c_str());
      tree->Branch((name+"PFMuon_GoodQuality").c_str() ,PFMuon_GoodQuality  ,(name+"PFMuon_GoodQuality["+name+"nPFMuon]/I").c_str());

      //--------------------------------------
      // secondary vertex information
      //--------------------------------------
      tree->Branch((name+"nSV").c_str()                ,&nSV               ,(name+"nSV/I").c_str());
      tree->Branch((name+"SV_x").c_str()               ,SV_x                 ,(name+"SV_x["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_y").c_str()               ,SV_y                 ,(name+"SV_y["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_z").c_str()               ,SV_z                 ,(name+"SV_z["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_ex").c_str()              ,SV_ex              ,(name+"SV_ex["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_ey").c_str()              ,SV_ey              ,(name+"SV_ey["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_ez").c_str()              ,SV_ez              ,(name+"SV_ez["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_chi2").c_str()            ,SV_chi2            ,(name+"SV_chi2["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_ndf").c_str()             ,SV_ndf             ,(name+"SV_ndf["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_flight").c_str()          ,SV_flight          ,(name+"SV_flight["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_flightErr").c_str()       ,SV_flightErr       ,(name+"SV_flightErr["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_deltaR_jet").c_str()      ,SV_deltaR_jet      ,(name+"SV_deltaR_jet["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_deltaR_sum_jet").c_str()  ,SV_deltaR_sum_jet  ,(name+"SV_deltaR_sum_jet["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_deltaR_sum_dir").c_str()  ,SV_deltaR_sum_dir  ,(name+"SV_deltaR_sum_dir["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_energy_ratio").c_str()    ,SV_energy_ratio    ,(name+"SV_energy_ratio["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_aboveC").c_str()          ,SV_aboveC             ,(name+"SV_aboveC["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_vtx_pt").c_str()          ,SV_vtx_pt             ,(name+"SV_vtx_pt["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_flight2D").c_str()        ,SV_flight2D        ,(name+"SV_flight2D["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_flight2DErr").c_str()     ,SV_flight2DErr     ,(name+"SV_flight2DErr["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_totCharge").c_str()       ,SV_totCharge       ,(name+"SV_totCharge ["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_vtxDistJetAxis").c_str()  ,SV_vtxDistJetAxis  ,(name+"SV_vtxDistJetAxis ["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_nTrk").c_str()            ,SV_nTrk            ,(name+"SV_nTrk["+name+"nSV]/I").c_str());
      tree->Branch((name+"SV_nTrk_firstVxt").c_str()   ,SV_nTrk_firstVxt   ,(name+"SV_nTrk_firstVxt["+name+"nSV]/I").c_str());
      tree->Branch((name+"SV_mass").c_str()            ,SV_mass            ,(name+"SV_mass["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_vtx_eta").c_str()         ,SV_vtx_eta         ,(name+"SV_vtx_eta["+name+"nSV]/F").c_str());
      tree->Branch((name+"SV_vtx_phi").c_str()         ,SV_vtx_phi         ,(name+"SV_vtx_phi["+name+"nSV]/F").c_str());
    }

    void RegisterJetTrackTree(TTree *tree, std::string name="") {
      if(name!="") name += ".";
      //--------------------------------------
      // track information
      //--------------------------------------
      tree->Branch((name+"nTrack").c_str()           ,&nTrack          ,(name+"nTrack/I").c_str());
      tree->Branch((name+"Track_dxy").c_str()        ,Track_dxy             ,(name+"Track_dxy["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_dz").c_str()         ,Track_dz         ,(name+"Track_dz["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_zIP").c_str()        ,Track_zIP             ,(name+"Track_zIP["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_length").c_str()     ,Track_length     ,(name+"Track_length["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_dist").c_str()       ,Track_dist            ,(name+"Track_dist["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_IP2D").c_str()       ,Track_IP2D            ,(name+"Track_IP2D["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_IP2Dsig").c_str()    ,Track_IP2Dsig    ,(name+"Track_IP2Dsig["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_IP2Derr").c_str()    ,Track_IP2Derr    ,(name+"Track_IP2Derr["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_IP").c_str()         ,Track_IP         ,(name+"Track_IP["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_IPsig").c_str()      ,Track_IPsig      ,(name+"Track_IPsig["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_IPerr").c_str()      ,Track_IPerr      ,(name+"Track_IPerr["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_Proba").c_str()      ,Track_Proba      ,(name+"Track_Proba["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_p").c_str()          ,Track_p          ,(name+"Track_p["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_pt").c_str()         ,Track_pt         ,(name+"Track_pt["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_eta").c_str()        ,Track_eta             ,(name+"Track_eta["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_phi").c_str()        ,Track_phi             ,(name+"Track_phi["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_chi2").c_str()       ,Track_chi2            ,(name+"Track_chi2["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_charge").c_str()     ,Track_charge     ,(name+"Track_charge["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_history").c_str()    ,Track_history    ,(name+"Track_history["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitStrip").c_str()  ,Track_nHitStrip  ,(name+"Track_nHitStrip["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitPixel").c_str()  ,Track_nHitPixel  ,(name+"Track_nHitPixel["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitAll").c_str()    ,Track_nHitAll    ,(name+"Track_nHitAll["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitTIB").c_str()    ,Track_nHitTIB    ,(name+"Track_nHitTIB["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitTID").c_str()    ,Track_nHitTID    ,(name+"Track_nHitTID["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitTOB").c_str()    ,Track_nHitTOB    ,(name+"Track_nHitTOB["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitTEC").c_str()    ,Track_nHitTEC    ,(name+"Track_nHitTEC["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitPXB").c_str()    ,Track_nHitPXB    ,(name+"Track_nHitPXB["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_nHitPXF").c_str()    ,Track_nHitPXF    ,(name+"Track_nHitPXF["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_isHitL1").c_str()    ,Track_isHitL1    ,(name+"Track_isHitL1["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_PV").c_str()         ,Track_PV         ,(name+"Track_PV["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_SV").c_str()         ,Track_SV         ,(name+"Track_SV["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_PVweight").c_str()   ,Track_PVweight   ,(name+"Track_PVweight["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_SVweight").c_str()   ,Track_SVweight   ,(name+"Track_SVweight["+name+"nTrack]/F").c_str());
      tree->Branch((name+"Track_isfromSV").c_str()   ,Track_isfromSV   ,(name+"Track_isfromSV["+name+"nTrack]/I").c_str());
      tree->Branch((name+"Track_category").c_str()   ,Track_category   ,(name+"Track_category["+name+"nTrack]/I").c_str());
    }

    void RegisterTagVarTree(TTree *tree, std::string name=""){
      if(name!="") name += ".";
      //--------------------------------------
      // TagInfo TaggingVariables
      //--------------------------------------
      tree->Branch((name+"Jet_nFirstTrkTagVar").c_str()  ,Jet_nFirstTrkTagVar  ,(name+"Jet_nFirstTrkTagVar["+name+"nJet]/I").c_str() );
      tree->Branch((name+"Jet_nLastTrkTagVar").c_str()   ,Jet_nLastTrkTagVar   ,(name+"Jet_nLastTrack["+name+"nJet]/I").c_str()      );
      tree->Branch((name+"Jet_nFirstSVTagVar").c_str()   ,Jet_nFirstSVTagVar   ,(name+"Jet_nFirstSVTagVar["+name+"nJet]/I").c_str()  );
      tree->Branch((name+"Jet_nLastSVTagVar").c_str()    ,Jet_nLastSVTagVar    ,(name+"Jet_nLastSVTagVar["+name+"nJet]/I").c_str()   );

      tree->Branch((name+"TagVar_jetNTracks").c_str()                   ,TagVar_jetNTracks                   ,(name+"TagVar_jetNTracks["+name+"nJet]/F").c_str()                   );
      tree->Branch((name+"TagVar_jetNSecondaryVertices").c_str()        ,TagVar_jetNSecondaryVertices        ,(name+"TagVar_jetNSecondaryVertices["+name+"nJet]/F").c_str()        );
      tree->Branch((name+"TagVar_chargedHadronEnergyFraction").c_str()  ,TagVar_chargedHadronEnergyFraction  ,(name+"TagVar_chargedHadronEnergyFraction["+name+"nJet]/F").c_str()  );
      tree->Branch((name+"TagVar_neutralHadronEnergyFraction").c_str()  ,TagVar_neutralHadronEnergyFraction  ,(name+"TagVar_neutralHadronEnergyFraction["+name+"nJet]/F").c_str()  );
      tree->Branch((name+"TagVar_photonEnergyFraction").c_str()         ,TagVar_photonEnergyFraction         ,(name+"TagVar_photonEnergyFraction["+name+"nJet]/F").c_str()         );
      tree->Branch((name+"TagVar_electronEnergyFraction").c_str()       ,TagVar_electronEnergyFraction       ,(name+"TagVar_electronEnergyFraction["+name+"nJet]/F").c_str()       );
      tree->Branch((name+"TagVar_muonEnergyFraction").c_str()           ,TagVar_muonEnergyFraction           ,(name+"TagVar_muonEnergyFraction["+name+"nJet]/F").c_str()           );
      tree->Branch((name+"TagVar_chargedHadronMultiplicity").c_str()    ,TagVar_chargedHadronMultiplicity    ,(name+"TagVar_chargedHadronMultiplicity["+name+"nJet]/F").c_str()    );
      tree->Branch((name+"TagVar_neutralHadronMultiplicity").c_str()    ,TagVar_neutralHadronMultiplicity    ,(name+"TagVar_neutralHadronMultiplicity["+name+"nJet]/F").c_str()    );
      tree->Branch((name+"TagVar_photonMultiplicity").c_str()           ,TagVar_photonMultiplicity           ,(name+"TagVar_photonMultiplicity["+name+"nJet]/F").c_str()           );
      tree->Branch((name+"TagVar_electronMultiplicity").c_str()         ,TagVar_electronMultiplicity         ,(name+"TagVar_electronMultiplicity["+name+"nJet]/F").c_str()         );
      tree->Branch((name+"TagVar_muonMultiplicity").c_str()             ,TagVar_muonMultiplicity             ,(name+"TagVar_muonMultiplicity["+name+"nJet]/F").c_str()             );

      tree->Branch((name+"nTrkTagVar").c_str()               ,&nTrkTagVar              ,(name+"nTrkTagVar/I").c_str()                                  );
      tree->Branch((name+"TagVar_trackMomentum").c_str()     ,TagVar_trackMomentum     ,(name+"TagVar_trackMomentum["+name+"nTrkTagVar]/F").c_str()    );
      tree->Branch((name+"TagVar_trackEta").c_str()          ,TagVar_trackEta          ,(name+"TagVar_trackEta["+name+"nTrkTagVar]/F").c_str()         );
      tree->Branch((name+"TagVar_trackPhi").c_str()          ,TagVar_trackPhi          ,(name+"TagVar_trackPhi["+name+"nTrkTagVar]/F").c_str()         );
      tree->Branch((name+"TagVar_trackPtRel").c_str()        ,TagVar_trackPtRel        ,(name+"TagVar_trackPtRel["+name+"nTrkTagVar]/F").c_str()       );
      tree->Branch((name+"TagVar_trackPPar").c_str()         ,TagVar_trackPPar         ,(name+"TagVar_trackPPar["+name+"nTrkTagVar]/F").c_str()        );
      tree->Branch((name+"TagVar_trackEtaRel").c_str()       ,TagVar_trackEtaRel       ,(name+"TagVar_trackEtaRel["+name+"nTrkTagVar]/F").c_str()      );
      tree->Branch((name+"TagVar_trackDeltaR").c_str()       ,TagVar_trackDeltaR       ,(name+"TagVar_trackDeltaR["+name+"nTrkTagVar]/F").c_str()      );
      tree->Branch((name+"TagVar_trackPtRatio").c_str()      ,TagVar_trackPtRatio      ,(name+"TagVar_trackPtRatio["+name+"nTrkTagVar]/F").c_str()     );
      tree->Branch((name+"TagVar_trackPParRatio").c_str()    ,TagVar_trackPParRatio    ,(name+"TagVar_trackPParRatio["+name+"nTrkTagVar]/F").c_str()   );
      tree->Branch((name+"TagVar_trackSip2dVal").c_str()     ,TagVar_trackSip2dVal     ,(name+"TagVar_trackSip2dVal["+name+"nTrkTagVar]/F").c_str()    );
      tree->Branch((name+"TagVar_trackSip2dSig").c_str()     ,TagVar_trackSip2dSig     ,(name+"TagVar_trackSip2dSig["+name+"nTrkTagVar]/F").c_str()    );
      tree->Branch((name+"TagVar_trackSip3dVal").c_str()     ,TagVar_trackSip3dVal     ,(name+"TagVar_trackSip3dVal["+name+"nTrkTagVar]/F").c_str()    );
      tree->Branch((name+"TagVar_trackSip3dSig").c_str()     ,TagVar_trackSip3dSig     ,(name+"TagVar_trackSip3dSig["+name+"nTrkTagVar]/F").c_str()    );
      tree->Branch((name+"TagVar_trackDecayLenVal").c_str()  ,TagVar_trackDecayLenVal  ,(name+"TagVar_trackDecayLenVal["+name+"nTrkTagVar]/F").c_str() );
      tree->Branch((name+"TagVar_trackDecayLenSig").c_str()  ,TagVar_trackDecayLenSig  ,(name+"TagVar_trackDecayLenSig["+name+"nTrkTagVar]/F").c_str() );
      tree->Branch((name+"TagVar_trackJetDistVal").c_str()   ,TagVar_trackJetDistVal   ,(name+"TagVar_trackJetDistVal["+name+"nTrkTagVar]/F").c_str()  );
      tree->Branch((name+"TagVar_trackJetDistSig").c_str()   ,TagVar_trackJetDistSig   ,(name+"TagVar_trackJetDistSig["+name+"nTrkTagVar]/F").c_str()  );
      tree->Branch((name+"TagVar_trackChi2").c_str()         ,TagVar_trackChi2         ,(name+"TagVar_trackChi2["+name+"nTrkTagVar]/F").c_str()        );
      tree->Branch((name+"TagVar_trackNTotalHits").c_str()   ,TagVar_trackNTotalHits   ,(name+"TagVar_trackNTotalHits["+name+"nTrkTagVar]/F").c_str()  );
      tree->Branch((name+"TagVar_trackNPixelHits").c_str()   ,TagVar_trackNPixelHits   ,(name+"TagVar_trackNPixelHits["+name+"nTrkTagVar]/F").c_str()  );

      tree->Branch((name+"nSVTagVar").c_str()                       ,&nSVTagVar                      ,(name+"nSVTagVar/I").c_str()                                     );
      tree->Branch((name+"TagVar_vertexMass").c_str()               ,TagVar_vertexMass               ,(name+"TagVar_vertexMass["+name+"nSVTagVar]/F").c_str()          );
      tree->Branch((name+"TagVar_vertexNTracks").c_str()            ,TagVar_vertexNTracks            ,(name+"TagVar_vertexNTracks["+name+"nSVTagVar]/F").c_str()       );
      tree->Branch((name+"TagVar_vertexJetDeltaR").c_str()          ,TagVar_vertexJetDeltaR          ,(name+"TagVar_vertexJetDeltaR["+name+"nSVTagVar]/F").c_str()     );
      tree->Branch((name+"TagVar_flightDistance2dVal").c_str()      ,TagVar_flightDistance2dVal      ,(name+"TagVar_flightDistance2dVal["+name+"nSVTagVar]/F").c_str() );
      tree->Branch((name+"TagVar_flightDistance2dSig").c_str()      ,TagVar_flightDistance2dSig      ,(name+"TagVar_flightDistance2dSig["+name+"nSVTagVar]/F").c_str() );
      tree->Branch((name+"TagVar_flightDistance3dVal").c_str()      ,TagVar_flightDistance3dVal      ,(name+"TagVar_flightDistance3dVal["+name+"nSVTagVar]/F").c_str() );
      tree->Branch((name+"TagVar_flightDistance3dSig").c_str()      ,TagVar_flightDistance3dSig      ,(name+"TagVar_flightDistance3dSig["+name+"nSVTagVar]/F").c_str() );

    }

    void RegisterCSVTagVarTree(TTree *tree, std::string name=""){
      if(name!="") name += ".";
      //--------------------------------------
      // CSV TaggingVariables
      //--------------------------------------
      tree->Branch((name+"Jet_nFirstTrkTagVarCSV").c_str()        ,Jet_nFirstTrkTagVarCSV        ,(name+"Jet_nFirstTrkTagVarCSV["+name+"nJet]/I").c_str()       );
      tree->Branch((name+"Jet_nLastTrkTagVarCSV").c_str()         ,Jet_nLastTrkTagVarCSV         ,(name+"Jet_nLastTrackCSV["+name+"nJet]/I").c_str()            );
      tree->Branch((name+"Jet_nFirstTrkEtaRelTagVarCSV").c_str()  ,Jet_nFirstTrkEtaRelTagVarCSV  ,(name+"Jet_nFirstTrkEtaRelTagVarCSV["+name+"nJet]/I").c_str() );
      tree->Branch((name+"Jet_nLastTrkEtaRelTagVarCSV").c_str()   ,Jet_nLastTrkEtaRelTagVarCSV   ,(name+"Jet_nLastEtaRelTrackCSV["+name+"nJet]/I").c_str()      );

      tree->Branch((name+"TagVarCSV_trackJetPt").c_str()               ,TagVarCSV_trackJetPt               ,(name+"TagVarCSV_trackJetPt["+name+"nJet]/F").c_str()              );
      tree->Branch((name+"TagVarCSV_jetNTracks").c_str()               ,TagVarCSV_jetNTracks               ,(name+"TagVarCSV_jetNTracks["+name+"nJet]/F").c_str()              );
      tree->Branch((name+"TagVarCSV_jetNTracksEtaRel").c_str()         ,TagVarCSV_jetNTracksEtaRel         ,(name+"TagVarCSV_jetNTracksEtaRel["+name+"nJet]/F").c_str()        );
      tree->Branch((name+"TagVarCSV_trackSumJetEtRatio").c_str()       ,TagVarCSV_trackSumJetEtRatio       ,(name+"TagVarCSV_trackSumJetEtRatio["+name+"nJet]/F").c_str()      );
      tree->Branch((name+"TagVarCSV_trackSumJetDeltaR").c_str()        ,TagVarCSV_trackSumJetDeltaR        ,(name+"TagVarCSV_trackSumJetDeltaR["+name+"nJet]/F").c_str()       );
      tree->Branch((name+"TagVarCSV_trackSip2dValAboveCharm").c_str()  ,TagVarCSV_trackSip2dValAboveCharm  ,(name+"TagVarCSV_trackSip2dValAboveCharm["+name+"nJet]/F").c_str() );
      tree->Branch((name+"TagVarCSV_trackSip2dSigAboveCharm").c_str()  ,TagVarCSV_trackSip2dSigAboveCharm  ,(name+"TagVarCSV_trackSip2dSigAboveCharm["+name+"nJet]/F").c_str() );
      tree->Branch((name+"TagVarCSV_trackSip3dValAboveCharm").c_str()  ,TagVarCSV_trackSip3dValAboveCharm  ,(name+"TagVarCSV_trackSip3dValAboveCharm["+name+"nJet]/F").c_str() );
      tree->Branch((name+"TagVarCSV_trackSip3dSigAboveCharm").c_str()  ,TagVarCSV_trackSip3dSigAboveCharm  ,(name+"TagVarCSV_trackSip3dSigAboveCharm["+name+"nJet]/F").c_str() );
      tree->Branch((name+"TagVarCSV_vertexCategory").c_str()           ,TagVarCSV_vertexCategory           ,(name+"TagVarCSV_vertexCategory["+name+"nJet]/F").c_str()          );
      tree->Branch((name+"TagVarCSV_jetNSecondaryVertices").c_str()    ,TagVarCSV_jetNSecondaryVertices    ,(name+"TagVarCSV_jetNSecondaryVertices["+name+"nJet]/F").c_str()   );
      tree->Branch((name+"TagVarCSV_vertexMass").c_str()               ,TagVarCSV_vertexMass               ,(name+"TagVarCSV_vertexMass["+name+"nJet]/F").c_str()              );
      tree->Branch((name+"TagVarCSV_vertexNTracks").c_str()            ,TagVarCSV_vertexNTracks            ,(name+"TagVarCSV_vertexNTracks["+name+"nJet]/F").c_str()           );
      tree->Branch((name+"TagVarCSV_vertexEnergyRatio").c_str()        ,TagVarCSV_vertexEnergyRatio        ,(name+"TagVarCSV_vertexEnergyRatio["+name+"nJet]/F").c_str()       );
      tree->Branch((name+"TagVarCSV_vertexJetDeltaR").c_str()          ,TagVarCSV_vertexJetDeltaR          ,(name+"TagVarCSV_vertexJetDeltaR["+name+"nJet]/F").c_str()         );
      tree->Branch((name+"TagVarCSV_flightDistance2dVal").c_str()      ,TagVarCSV_flightDistance2dVal      ,(name+"TagVarCSV_flightDistance2dVal["+name+"nJet]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_flightDistance2dSig").c_str()      ,TagVarCSV_flightDistance2dSig      ,(name+"TagVarCSV_flightDistance2dSig["+name+"nJet]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_flightDistance3dVal").c_str()      ,TagVarCSV_flightDistance3dVal      ,(name+"TagVarCSV_flightDistance3dVal["+name+"nJet]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_flightDistance3dSig").c_str()      ,TagVarCSV_flightDistance3dSig      ,(name+"TagVarCSV_flightDistance3dSig["+name+"nJet]/F").c_str()     );

      tree->Branch((name+"nTrkTagVarCSV").c_str()               ,&nTrkTagVarCSV              ,(name+"nTrkTagVarCSV/I").c_str()                                      );
      tree->Branch((name+"nTrkEtaRelTagVarCSV").c_str()         ,&nTrkEtaRelTagVarCSV        ,(name+"nTrkEtaRelTagVarCSV/I").c_str()                                );
      tree->Branch((name+"TagVarCSV_trackMomentum").c_str()     ,TagVarCSV_trackMomentum     ,(name+"TagVarCSV_trackMomentum["+name+"nTrkTagVarCSV]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_trackEta").c_str()          ,TagVarCSV_trackEta          ,(name+"TagVarCSV_trackEta["+name+"nTrkTagVarCSV]/F").c_str()          );
      tree->Branch((name+"TagVarCSV_trackPhi").c_str()          ,TagVarCSV_trackPhi          ,(name+"TagVarCSV_trackPhi["+name+"nTrkTagVarCSV]/F").c_str()          );
      tree->Branch((name+"TagVarCSV_trackPtRel").c_str()        ,TagVarCSV_trackPtRel        ,(name+"TagVarCSV_trackPtRel["+name+"nTrkTagVarCSV]/F").c_str()        );
      tree->Branch((name+"TagVarCSV_trackPPar").c_str()         ,TagVarCSV_trackPPar         ,(name+"TagVarCSV_trackPPar["+name+"nTrkTagVarCSV]/F").c_str()         );
      tree->Branch((name+"TagVarCSV_trackDeltaR").c_str()       ,TagVarCSV_trackDeltaR       ,(name+"TagVarCSV_trackDeltaR["+name+"nTrkTagVarCSV]/F").c_str()       );
      tree->Branch((name+"TagVarCSV_trackPtRatio").c_str()      ,TagVarCSV_trackPtRatio      ,(name+"TagVarCSV_trackPtRatio["+name+"nTrkTagVarCSV]/F").c_str()      );
      tree->Branch((name+"TagVarCSV_trackPParRatio").c_str()    ,TagVarCSV_trackPParRatio    ,(name+"TagVarCSV_trackPParRatio["+name+"nTrkTagVarCSV]/F").c_str()    );
      tree->Branch((name+"TagVarCSV_trackSip2dVal").c_str()     ,TagVarCSV_trackSip2dVal     ,(name+"TagVarCSV_trackSip2dVal["+name+"nTrkTagVarCSV]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_trackSip2dSig").c_str()     ,TagVarCSV_trackSip2dSig     ,(name+"TagVarCSV_trackSip2dSig["+name+"nTrkTagVarCSV]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_trackSip3dVal").c_str()     ,TagVarCSV_trackSip3dVal     ,(name+"TagVarCSV_trackSip3dVal["+name+"nTrkTagVarCSV]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_trackSip3dSig").c_str()     ,TagVarCSV_trackSip3dSig     ,(name+"TagVarCSV_trackSip3dSig["+name+"nTrkTagVarCSV]/F").c_str()     );
      tree->Branch((name+"TagVarCSV_trackDecayLenVal").c_str()  ,TagVarCSV_trackDecayLenVal  ,(name+"TagVarCSV_trackDecayLenVal["+name+"nTrkTagVarCSV]/F").c_str()  );
      tree->Branch((name+"TagVarCSV_trackDecayLenSig").c_str()  ,TagVarCSV_trackDecayLenSig  ,(name+"TagVarCSV_trackDecayLenSig["+name+"nTrkTagVarCSV]/F").c_str()  );
      tree->Branch((name+"TagVarCSV_trackJetDistVal").c_str()   ,TagVarCSV_trackJetDistVal   ,(name+"TagVarCSV_trackJetDistVal["+name+"nTrkTagVarCSV]/F").c_str()   );
      tree->Branch((name+"TagVarCSV_trackJetDistSig").c_str()   ,TagVarCSV_trackJetDistSig   ,(name+"TagVarCSV_trackJetDistSig["+name+"nTrkTagVarCSV]/F").c_str()   );
      tree->Branch((name+"TagVarCSV_trackEtaRel").c_str()       ,TagVarCSV_trackEtaRel       ,(name+"TagVarCSV_trackEtaRel["+name+"nTrkEtaRelTagVarCSV]/F").c_str() );
    }

    void RegisterSubJetSpecificTree(TTree *tree, std::string name="") {
      if(name!="") name += ".";
      tree->Branch((name+"Jet_FatJetIdx").c_str(),   Jet_FatJetIdx   ,(name+"Jet_FatJetIdx["+name+"nJet]/I").c_str());
    }

    void RegisterFatJetSpecificTree(TTree *tree, std::string name="") {
      if(name!="") name += ".";
      tree->Branch((name+"Jet_ptGroomed").c_str(),   Jet_ptGroomed   ,(name+"Jet_ptGroomed["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_jesGroomed").c_str(),  Jet_jesGroomed  ,(name+"Jet_jesGroomed["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_etaGroomed").c_str(),  Jet_etaGroomed  ,(name+"Jet_etaGroomed["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_phiGroomed").c_str(),  Jet_phiGroomed  ,(name+"Jet_phiGroomed["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_massGroomed").c_str(), Jet_massGroomed ,(name+"Jet_massGroomed["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_tau1").c_str(),        Jet_tau1        ,(name+"Jet_tau1["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_tau2").c_str(),        Jet_tau2        ,(name+"Jet_tau2["+name+"nJet]/F").c_str());
      tree->Branch((name+"Jet_nSubJets").c_str(),    Jet_nSubJets    ,(name+"Jet_nSubJets["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nFirstSJ").c_str(),    Jet_nFirstSJ    ,(name+"Jet_nFirstSJ["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nLastSJ").c_str(),     Jet_nLastSJ     ,(name+"Jet_nLastSJ["+name+"nJet]/I").c_str());
      tree->Branch((name+"nSubJet").c_str(),         &nSubJet        ,(name+"nSubJet/I").c_str());
      tree->Branch((name+"SubJetIdx").c_str(),       SubJetIdx       ,(name+"SubJetIdx["+name+"nSubJet]/I").c_str());
      tree->Branch((name+"Jet_nsharedtracks").c_str(), Jet_nsharedtracks ,(name+"Jet_nsharedtracks["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nsubjettracks").c_str(), Jet_nsubjettracks ,(name+"Jet_nsubjettracks["+name+"nJet]/I").c_str());
      tree->Branch((name+"Jet_nsharedsubjettracks").c_str(), Jet_nsharedsubjettracks ,(name+"Jet_nsharedsubjettracks["+name+"nJet]/I").c_str());
    }

    //------------------------------------------------------------------------------------------------------------------

    void ReadTree(TTree *tree, std::string name="") {
      if (name!="") name += ".";
      tree->SetBranchAddress((name+"nJet").c_str(),            &nJet           );
      tree->SetBranchAddress((name+"Jet_pt").c_str(),          Jet_pt          );
      tree->SetBranchAddress((name+"Jet_genpt").c_str(),       Jet_genpt       );
      tree->SetBranchAddress((name+"Jet_residual").c_str(),    Jet_residual    );
      tree->SetBranchAddress((name+"Jet_jes").c_str(),         Jet_jes         );
      tree->SetBranchAddress((name+"Jet_eta").c_str(),         Jet_eta         );
      tree->SetBranchAddress((name+"Jet_phi").c_str(),         Jet_phi         );
      tree->SetBranchAddress((name+"Jet_mass").c_str(),        Jet_mass        );
      tree->SetBranchAddress((name+"Jet_ntracks").c_str(),     Jet_ntracks     );
      tree->SetBranchAddress((name+"Jet_nseltracks").c_str(),  Jet_nseltracks  );
      tree->SetBranchAddress((name+"Jet_flavour").c_str(),     Jet_flavour     );
      tree->SetBranchAddress((name+"Jet_Ip2N").c_str(),        Jet_Ip2N        );
      tree->SetBranchAddress((name+"Jet_Ip2P").c_str(),        Jet_Ip2P        );
      tree->SetBranchAddress((name+"Jet_Ip3N").c_str(),        Jet_Ip3N        );
      tree->SetBranchAddress((name+"Jet_Ip3P").c_str(),        Jet_Ip3P        );
      tree->SetBranchAddress((name+"Jet_ProbaN").c_str(),      Jet_ProbaN      );
      tree->SetBranchAddress((name+"Jet_ProbaP").c_str(),      Jet_ProbaP      );
      tree->SetBranchAddress((name+"Jet_Proba").c_str(),       Jet_Proba       );
      tree->SetBranchAddress((name+"Jet_BprobN").c_str(),      Jet_BprobN      );
      tree->SetBranchAddress((name+"Jet_BprobP").c_str(),      Jet_BprobP      );
      tree->SetBranchAddress((name+"Jet_Bprob").c_str(),       Jet_Bprob       );
      tree->SetBranchAddress((name+"Jet_SvxN").c_str(),        Jet_SvxN        );
      tree->SetBranchAddress((name+"Jet_Svx").c_str(),         Jet_Svx         );
      tree->SetBranchAddress((name+"Jet_SvxNHP").c_str(),      Jet_SvxNHP      );
      tree->SetBranchAddress((name+"Jet_SvxHP").c_str(),       Jet_SvxHP       );
      tree->SetBranchAddress((name+"Jet_SvxMass").c_str(),     Jet_SvxMass     );
      tree->SetBranchAddress((name+"Jet_CombSvxN").c_str(),    Jet_CombSvxN    );
      tree->SetBranchAddress((name+"Jet_CombSvxP").c_str(),    Jet_CombSvxP    );
      tree->SetBranchAddress((name+"Jet_CombSvx").c_str(),     Jet_CombSvx     );

      tree->SetBranchAddress((name+"Jet_RetCombSvxN").c_str(), Jet_RetCombSvxN    );
      tree->SetBranchAddress((name+"Jet_RetCombSvxP").c_str(), Jet_RetCombSvxP    );
      tree->SetBranchAddress((name+"Jet_RetCombSvx").c_str(), Jet_RetCombSvx     );

      tree->SetBranchAddress((name+"Jet_CombCSVJP_N").c_str(), Jet_CombCSVJP_N    );
      tree->SetBranchAddress((name+"Jet_CombCSVJP_P").c_str(), Jet_CombCSVJP_P    );
      tree->SetBranchAddress((name+"Jet_CombCSVJP").c_str(),    Jet_CombCSVJP      );

      tree->SetBranchAddress((name+"Jet_CombCSVSL_N").c_str(), Jet_CombCSVSL_N    );
      tree->SetBranchAddress((name+"Jet_CombCSVSL_P").c_str(), Jet_CombCSVSL_P    );
      tree->SetBranchAddress((name+"Jet_CombCSVSL").c_str(), Jet_CombCSVSL      );

      tree->SetBranchAddress((name+"Jet_CombCSVJPSL_N").c_str(), Jet_CombCSVJPSL_N  );
      tree->SetBranchAddress((name+"Jet_CombCSVJPSL_P").c_str(), Jet_CombCSVJPSL_P  );
      tree->SetBranchAddress((name+"Jet_CombCSVJPSL").c_str(), Jet_CombCSVJPSL    );

      // tree->SetBranchAddress((name+"Jet_SimpIVF_HP").c_str(),  Jet_SimpIVF_HP  );
      // tree->SetBranchAddress((name+"Jet_SimpIVF_HE").c_str(),  Jet_SimpIVF_HE  );
      // tree->SetBranchAddress((name+"Jet_DoubIVF_HE").c_str(),  Jet_DoubIVF_HE  );
      // tree->SetBranchAddress((name+"Jet_CombIVF").c_str(),     Jet_CombIVF     );
      // tree->SetBranchAddress((name+"Jet_CombIVF_P").c_str(), Jet_CombIVF_P   );

      tree->SetBranchAddress((name+"Jet_SoftMuN").c_str(), Jet_SoftMuN     );
      tree->SetBranchAddress((name+"Jet_SoftMuP").c_str(), Jet_SoftMuP     );
      tree->SetBranchAddress((name+"Jet_SoftMu").c_str(),  Jet_SoftMu      );

      tree->SetBranchAddress((name+"Jet_SoftElN").c_str(),     Jet_SoftElN     );
      tree->SetBranchAddress((name+"Jet_SoftElP").c_str(),     Jet_SoftElP     );
      tree->SetBranchAddress((name+"Jet_SoftEl").c_str(),      Jet_SoftEl      );

      tree->SetBranchAddress((name+"Jet_hist1").c_str(),       Jet_hist1       );
      tree->SetBranchAddress((name+"Jet_hist2").c_str(),       Jet_hist2       );
      tree->SetBranchAddress((name+"Jet_hist3").c_str(),       Jet_hist3       );
      tree->SetBranchAddress((name+"Jet_histJet").c_str(),     Jet_histJet     );
      tree->SetBranchAddress((name+"Jet_histSvx").c_str(),     Jet_histSvx     );

      tree->SetBranchAddress((name+"Jet_nFirstTrack").c_str(), Jet_nFirstTrack );
      tree->SetBranchAddress((name+"Jet_nLastTrack").c_str(),  Jet_nLastTrack  );
      tree->SetBranchAddress((name+"Jet_nFirstSV").c_str(),    Jet_nFirstSV    );
      tree->SetBranchAddress((name+"Jet_nLastSV").c_str(),     Jet_nLastSV     );
      tree->SetBranchAddress((name+"Jet_SV_multi").c_str(),    Jet_SV_multi      );
      tree->SetBranchAddress((name+"Jet_nFirstTrkInc").c_str(),Jet_nFirstTrkInc );
      tree->SetBranchAddress((name+"Jet_nLastTrkInc").c_str(), Jet_nLastTrkInc  );

      tree->SetBranchAddress((name+"Jet_VtxCat").c_str(),      Jet_VtxCat  );
      tree->SetBranchAddress((name+"Jet_looseID").c_str(),     Jet_looseID);
      tree->SetBranchAddress((name+"Jet_tightID").c_str(),     Jet_tightID);

      //--------------------------------------
      // Inclusive Track information for PtRel template
      //--------------------------------------
      tree->SetBranchAddress((name+"nTrkInc").c_str()      ,&nTrkInc     ) ;
      tree->SetBranchAddress((name+"TrkInc_pt").c_str()    ,TrkInc_pt    ) ;
      tree->SetBranchAddress((name+"TrkInc_eta").c_str()   ,TrkInc_eta   ) ;
      tree->SetBranchAddress((name+"TrkInc_phi").c_str()   ,TrkInc_phi   ) ;
      tree->SetBranchAddress((name+"TrkInc_ptrel").c_str() ,TrkInc_ptrel ) ;
      tree->SetBranchAddress((name+"TrkInc_IPsig").c_str() ,TrkInc_IPsig ) ;
      tree->SetBranchAddress((name+"TrkInc_IP").c_str()    ,TrkInc_IP    ) ;

      //--------------------------------------
      // muon information
      //--------------------------------------
      tree->SetBranchAddress((name+"nMuon").c_str()         ,&nMuon        ) ;
      tree->SetBranchAddress((name+"Muon_IdxJet").c_str()   ,Muon_IdxJet   ) ;
      tree->SetBranchAddress((name+"Muon_nMuHit").c_str()   ,Muon_nMuHit   ) ;
      tree->SetBranchAddress((name+"Muon_nTkHit").c_str()   ,Muon_nTkHit   ) ;
      tree->SetBranchAddress((name+"Muon_nPixHit").c_str()  ,Muon_nPixHit  ) ;
      tree->SetBranchAddress((name+"Muon_nOutHit").c_str()  ,Muon_nOutHit  ) ;
      tree->SetBranchAddress((name+"Muon_isGlobal").c_str() ,Muon_isGlobal ) ;
      tree->SetBranchAddress((name+"Muon_nMatched").c_str() ,Muon_nMatched ) ;
      tree->SetBranchAddress((name+"Muon_chi2").c_str()     ,Muon_chi2     ) ;
      tree->SetBranchAddress((name+"Muon_chi2Tk").c_str()   ,Muon_chi2Tk   ) ;
      tree->SetBranchAddress((name+"Muon_pt").c_str()       ,Muon_pt       ) ;
      tree->SetBranchAddress((name+"Muon_eta").c_str()      ,Muon_eta      ) ;
      tree->SetBranchAddress((name+"Muon_phi").c_str()      ,Muon_phi      ) ;
      tree->SetBranchAddress((name+"Muon_ptrel").c_str()    ,Muon_ptrel    ) ;
      tree->SetBranchAddress((name+"Muon_vz").c_str()	    ,Muon_vz	   ) ;
      tree->SetBranchAddress((name+"Muon_hist").c_str()     ,Muon_hist     ) ;
      tree->SetBranchAddress((name+"Muon_TrackIdx").c_str() ,Muon_TrackIdx ) ;
      tree->SetBranchAddress((name+"Muon_IPsig").c_str()    ,Muon_IPsig    ) ;
      tree->SetBranchAddress((name+"Muon_IP").c_str()       ,Muon_IP       ) ;
      tree->SetBranchAddress((name+"Muon_IP2Dsig").c_str()  ,Muon_IP2Dsig  ) ;
      tree->SetBranchAddress((name+"Muon_IP2D").c_str()     ,Muon_IP2D     ) ;
      tree->SetBranchAddress((name+"Muon_Proba").c_str()    ,Muon_Proba    ) ;
      tree->SetBranchAddress((name+"Muon_deltaR").c_str()   ,Muon_deltaR   ) ;
      tree->SetBranchAddress((name+"Muon_ratio").c_str()    ,Muon_ratio    ) ;
      tree->SetBranchAddress((name+"Muon_ratioRel").c_str() ,Muon_ratioRel ) ;

      //--------------------------------------
      // pf electron information
      //--------------------------------------
      tree->SetBranchAddress((name+"nPFElectron").c_str()         ,&nPFElectron        ) ;
      tree->SetBranchAddress((name+"PFElectron_IdxJet").c_str()   ,PFElectron_IdxJet  ) ;
      tree->SetBranchAddress((name+"PFElectron_pt").c_str()       ,PFElectron_pt      ) ;
      tree->SetBranchAddress((name+"PFElectron_eta").c_str()      ,PFElectron_eta     ) ;
      tree->SetBranchAddress((name+"PFElectron_phi").c_str()      ,PFElectron_phi     ) ;
      tree->SetBranchAddress((name+"PFElectron_ptrel").c_str()    ,PFElectron_ptrel   ) ;
      tree->SetBranchAddress((name+"PFElectron_deltaR").c_str()   ,PFElectron_deltaR  ) ;
      tree->SetBranchAddress((name+"PFElectron_ratio").c_str()    ,PFElectron_ratio   ) ;
      tree->SetBranchAddress((name+"PFElectron_ratioRel").c_str() ,PFElectron_ratioRel) ;
      tree->SetBranchAddress((name+"PFElectron_IP").c_str()       ,PFElectron_IP      ) ;
      tree->SetBranchAddress((name+"PFElectron_IP2D").c_str()     ,PFElectron_IP2D    ) ;

      //--------------------------------------
      // pf muon information
      //--------------------------------------
      tree->SetBranchAddress((name+"nPFMuon").c_str()            ,&nPFMuon            ) ;
      tree->SetBranchAddress((name+"PFMuon_IdxJet").c_str()      ,PFMuon_IdxJet       ) ;
      tree->SetBranchAddress((name+"PFMuon_pt").c_str()          ,PFMuon_pt           ) ;
      tree->SetBranchAddress((name+"PFMuon_eta").c_str()         ,PFMuon_eta          ) ;
      tree->SetBranchAddress((name+"PFMuon_phi").c_str()         ,PFMuon_phi          ) ;
      tree->SetBranchAddress((name+"PFMuon_ptrel").c_str()       ,PFMuon_ptrel        ) ;
      tree->SetBranchAddress((name+"PFMuon_deltaR").c_str()      ,PFMuon_deltaR       ) ;
      tree->SetBranchAddress((name+"PFMuon_ratio").c_str()       ,PFMuon_ratio        ) ;
      tree->SetBranchAddress((name+"PFMuon_ratioRel").c_str()    ,PFMuon_ratioRel     ) ;
      tree->SetBranchAddress((name+"PFMuon_IP").c_str()          ,PFMuon_IP           ) ;
      tree->SetBranchAddress((name+"PFMuon_IP2D").c_str()        ,PFMuon_IP2D         ) ;
      tree->SetBranchAddress((name+"PFMuon_GoodQuality").c_str() ,PFMuon_GoodQuality  ) ;

      //--------------------------------------
      // secondary vertex information
      //--------------------------------------
      tree->SetBranchAddress((name+"nSV").c_str()                ,&nSV               ) ;
      tree->SetBranchAddress((name+"SV_x").c_str()               ,SV_x                     ) ;
      tree->SetBranchAddress((name+"SV_y").c_str()               ,SV_y                     ) ;
      tree->SetBranchAddress((name+"SV_z").c_str()               ,SV_z                     ) ;
      tree->SetBranchAddress((name+"SV_ex").c_str()              ,SV_ex                  ) ;
      tree->SetBranchAddress((name+"SV_ey").c_str()              ,SV_ey                  ) ;
      tree->SetBranchAddress((name+"SV_ez").c_str()              ,SV_ez                  ) ;
      tree->SetBranchAddress((name+"SV_chi2").c_str()            ,SV_chi2            ) ;
      tree->SetBranchAddress((name+"SV_ndf").c_str()             ,SV_ndf                 ) ;
      tree->SetBranchAddress((name+"SV_flight").c_str()          ,SV_flight          ) ;
      tree->SetBranchAddress((name+"SV_flightErr").c_str()       ,SV_flightErr       ) ;
      tree->SetBranchAddress((name+"SV_deltaR_jet").c_str()      ,SV_deltaR_jet      ) ;
      tree->SetBranchAddress((name+"SV_deltaR_sum_jet").c_str()  ,SV_deltaR_sum_jet  ) ;
      tree->SetBranchAddress((name+"SV_deltaR_sum_dir").c_str()  ,SV_deltaR_sum_dir  ) ;
      tree->SetBranchAddress((name+"SV_energy_ratio").c_str()    ,SV_energy_ratio    ) ;
      tree->SetBranchAddress((name+"SV_aboveC").c_str()          ,SV_aboveC          ) ;
      tree->SetBranchAddress((name+"SV_vtx_pt").c_str()          ,SV_vtx_pt          ) ;
      tree->SetBranchAddress((name+"SV_flight2D").c_str()        ,SV_flight2D        ) ;
      tree->SetBranchAddress((name+"SV_flight2DErr").c_str()     ,SV_flight2DErr     ) ;
      tree->SetBranchAddress((name+"SV_totCharge").c_str()       ,SV_totCharge       ) ;
      tree->SetBranchAddress((name+"SV_vtxDistJetAxis").c_str()  ,SV_vtxDistJetAxis  ) ;
      tree->SetBranchAddress((name+"SV_nTrk").c_str()            ,SV_nTrk            ) ;
      tree->SetBranchAddress((name+"SV_nTrk_firstVxt").c_str()   ,SV_nTrk_firstVxt   ) ;
      tree->SetBranchAddress((name+"SV_mass").c_str()            ,SV_mass            ) ;
      tree->SetBranchAddress((name+"SV_vtx_eta").c_str()         ,SV_vtx_eta         ) ;
      tree->SetBranchAddress((name+"SV_vtx_phi").c_str()         ,SV_vtx_phi         ) ;
    }

    void ReadJetTrackTree(TTree *tree, std::string name="") {
      if (name!="") name += ".";
      //--------------------------------------
      // track information
      //--------------------------------------
      tree->SetBranchAddress((name+"nTrack").c_str()          ,&nTrack            ) ;
      tree->SetBranchAddress((name+"Track_dxy").c_str()       ,Track_dxy          ) ;
      tree->SetBranchAddress((name+"Track_dz").c_str()        ,Track_dz           ) ;
      tree->SetBranchAddress((name+"Track_zIP").c_str()       ,Track_zIP          ) ;
      tree->SetBranchAddress((name+"Track_length").c_str()    ,Track_length   ) ;
      tree->SetBranchAddress((name+"Track_dist").c_str()      ,Track_dist     ) ;
      tree->SetBranchAddress((name+"Track_IP2D").c_str()      ,Track_IP2D     ) ;
      tree->SetBranchAddress((name+"Track_IP2Dsig").c_str()   ,Track_IP2Dsig  ) ;
      tree->SetBranchAddress((name+"Track_IP2Derr").c_str()   ,Track_IP2Derr  ) ;
      tree->SetBranchAddress((name+"Track_IP").c_str()        ,Track_IP           ) ;
      tree->SetBranchAddress((name+"Track_IPsig").c_str()     ,Track_IPsig    ) ;
      tree->SetBranchAddress((name+"Track_IPerr").c_str()     ,Track_IPerr    ) ;
      tree->SetBranchAddress((name+"Track_Proba").c_str()     ,Track_Proba    ) ;
      tree->SetBranchAddress((name+"Track_p").c_str()         ,Track_p            ) ;
      tree->SetBranchAddress((name+"Track_pt").c_str()        ,Track_pt           ) ;
      tree->SetBranchAddress((name+"Track_eta").c_str()       ,Track_eta          ) ;
      tree->SetBranchAddress((name+"Track_phi").c_str()       ,Track_phi          ) ;
      tree->SetBranchAddress((name+"Track_chi2").c_str()      ,Track_chi2     ) ;
      tree->SetBranchAddress((name+"Track_charge").c_str()    ,Track_charge   ) ;
      tree->SetBranchAddress((name+"Track_history").c_str()   ,Track_history  ) ;
      tree->SetBranchAddress((name+"Track_nHitStrip").c_str() ,Track_nHitStrip) ;
      tree->SetBranchAddress((name+"Track_nHitPixel").c_str() ,Track_nHitPixel) ;
      tree->SetBranchAddress((name+"Track_nHitAll").c_str()   ,Track_nHitAll  ) ;
      tree->SetBranchAddress((name+"Track_nHitTIB").c_str()   ,Track_nHitTIB  ) ;
      tree->SetBranchAddress((name+"Track_nHitTID").c_str()   ,Track_nHitTID  ) ;
      tree->SetBranchAddress((name+"Track_nHitTOB").c_str()   ,Track_nHitTOB  ) ;
      tree->SetBranchAddress((name+"Track_nHitTEC").c_str()   ,Track_nHitTEC  ) ;
      tree->SetBranchAddress((name+"Track_nHitPXB").c_str()   ,Track_nHitPXB  ) ;
      tree->SetBranchAddress((name+"Track_nHitPXF").c_str()   ,Track_nHitPXF  ) ;
      tree->SetBranchAddress((name+"Track_isHitL1").c_str()   ,Track_isHitL1  ) ;
      tree->SetBranchAddress((name+"Track_PV").c_str()        ,Track_PV       ) ;
      tree->SetBranchAddress((name+"Track_SV").c_str()        ,Track_SV       ) ;
      tree->SetBranchAddress((name+"Track_PVweight").c_str()  ,Track_PVweight ) ;
      tree->SetBranchAddress((name+"Track_SVweight").c_str()  ,Track_SVweight ) ;
      tree->SetBranchAddress((name+"Track_isfromSV").c_str()  ,Track_isfromSV ) ;
      tree->SetBranchAddress((name+"Track_category").c_str()  ,Track_category ) ;
    }

    void ReadTagVarTree(TTree *tree, std::string name=""){
      if (name!="") name += ".";
      //--------------------------------------
      // TagInfo TaggingVariables
      //--------------------------------------
      tree->SetBranchAddress((name+"Jet_nFirstTrkTagVar").c_str()  ,Jet_nFirstTrkTagVar );
      tree->SetBranchAddress((name+"Jet_nLastTrkTagVar").c_str()   ,Jet_nLastTrkTagVar  );
      tree->SetBranchAddress((name+"Jet_nFirstSVTagVar").c_str()   ,Jet_nFirstSVTagVar  );
      tree->SetBranchAddress((name+"Jet_nLastSVTagVar").c_str()    ,Jet_nLastSVTagVar   );

      tree->SetBranchAddress((name+"TagVar_jetNTracks").c_str()                   ,TagVar_jetNTracks                  );
      tree->SetBranchAddress((name+"TagVar_jetNSecondaryVertices").c_str()        ,TagVar_jetNSecondaryVertices       );
      tree->SetBranchAddress((name+"TagVar_chargedHadronEnergyFraction").c_str()  ,TagVar_chargedHadronEnergyFraction );
      tree->SetBranchAddress((name+"TagVar_neutralHadronEnergyFraction").c_str()  ,TagVar_neutralHadronEnergyFraction );
      tree->SetBranchAddress((name+"TagVar_photonEnergyFraction").c_str()         ,TagVar_photonEnergyFraction        );
      tree->SetBranchAddress((name+"TagVar_electronEnergyFraction").c_str()       ,TagVar_electronEnergyFraction      );
      tree->SetBranchAddress((name+"TagVar_muonEnergyFraction").c_str()           ,TagVar_muonEnergyFraction          );
      tree->SetBranchAddress((name+"TagVar_chargedHadronMultiplicity").c_str()    ,TagVar_chargedHadronMultiplicity   );
      tree->SetBranchAddress((name+"TagVar_neutralHadronMultiplicity").c_str()    ,TagVar_neutralHadronMultiplicity   );
      tree->SetBranchAddress((name+"TagVar_photonMultiplicity").c_str()           ,TagVar_photonMultiplicity          );
      tree->SetBranchAddress((name+"TagVar_electronMultiplicity").c_str()         ,TagVar_electronMultiplicity        );
      tree->SetBranchAddress((name+"TagVar_muonMultiplicity").c_str()             ,TagVar_muonMultiplicity            );

      tree->SetBranchAddress((name+"nTrkTagVar").c_str()               ,&nTrkTagVar             );
      tree->SetBranchAddress((name+"TagVar_trackMomentum").c_str()     ,TagVar_trackMomentum    );
      tree->SetBranchAddress((name+"TagVar_trackEta").c_str()          ,TagVar_trackEta         );
      tree->SetBranchAddress((name+"TagVar_trackPhi").c_str()          ,TagVar_trackPhi         );
      tree->SetBranchAddress((name+"TagVar_trackPtRel").c_str()        ,TagVar_trackPtRel       );
      tree->SetBranchAddress((name+"TagVar_trackPPar").c_str()         ,TagVar_trackPPar        );
      tree->SetBranchAddress((name+"TagVar_trackEtaRel").c_str()       ,TagVar_trackEtaRel      );
      tree->SetBranchAddress((name+"TagVar_trackDeltaR").c_str()       ,TagVar_trackDeltaR      );
      tree->SetBranchAddress((name+"TagVar_trackPtRatio").c_str()      ,TagVar_trackPtRatio     );
      tree->SetBranchAddress((name+"TagVar_trackPParRatio").c_str()    ,TagVar_trackPParRatio   );
      tree->SetBranchAddress((name+"TagVar_trackSip2dVal").c_str()     ,TagVar_trackSip2dVal    );
      tree->SetBranchAddress((name+"TagVar_trackSip2dSig").c_str()     ,TagVar_trackSip2dSig    );
      tree->SetBranchAddress((name+"TagVar_trackSip3dVal").c_str()     ,TagVar_trackSip3dVal    );
      tree->SetBranchAddress((name+"TagVar_trackSip3dSig").c_str()     ,TagVar_trackSip3dSig    );
      tree->SetBranchAddress((name+"TagVar_trackDecayLenVal").c_str()  ,TagVar_trackDecayLenVal );
      tree->SetBranchAddress((name+"TagVar_trackDecayLenSig").c_str()  ,TagVar_trackDecayLenSig );
      tree->SetBranchAddress((name+"TagVar_trackJetDistVal").c_str()   ,TagVar_trackJetDistVal  );
      tree->SetBranchAddress((name+"TagVar_trackJetDistSig").c_str()   ,TagVar_trackJetDistSig  );
      tree->SetBranchAddress((name+"TagVar_trackChi2").c_str()         ,TagVar_trackChi2        );
      tree->SetBranchAddress((name+"TagVar_trackNTotalHits").c_str()   ,TagVar_trackNTotalHits  );
      tree->SetBranchAddress((name+"TagVar_trackNPixelHits").c_str()   ,TagVar_trackNPixelHits  );

      tree->SetBranchAddress((name+"nSVTagVar").c_str()                       ,&nSVTagVar                     );
      tree->SetBranchAddress((name+"TagVar_vertexMass").c_str()               ,TagVar_vertexMass              );
      tree->SetBranchAddress((name+"TagVar_vertexNTracks").c_str()            ,TagVar_vertexNTracks           );
      tree->SetBranchAddress((name+"TagVar_vertexJetDeltaR").c_str()          ,TagVar_vertexJetDeltaR         );
      tree->SetBranchAddress((name+"TagVar_flightDistance2dVal").c_str()      ,TagVar_flightDistance2dVal     );
      tree->SetBranchAddress((name+"TagVar_flightDistance2dSig").c_str()      ,TagVar_flightDistance2dSig     );
      tree->SetBranchAddress((name+"TagVar_flightDistance3dVal").c_str()      ,TagVar_flightDistance3dVal     );
      tree->SetBranchAddress((name+"TagVar_flightDistance3dSig").c_str()      ,TagVar_flightDistance3dSig     );

    }

    void ReadCSVTagVarTree(TTree *tree, std::string name=""){
      if (name!="") name += ".";
      //--------------------------------------
      // CSV TaggingVariables
      //--------------------------------------
      tree->SetBranchAddress((name+"Jet_nFirstTrkTagVarCSV").c_str()        ,Jet_nFirstTrkTagVarCSV );
      tree->SetBranchAddress((name+"Jet_nLastTrkTagVarCSV").c_str()         ,Jet_nLastTrkTagVarCSV  );
      tree->SetBranchAddress((name+"Jet_nFirstTrkEtaRelTagVarCSV").c_str()  ,Jet_nFirstTrkEtaRelTagVarCSV );
      tree->SetBranchAddress((name+"Jet_nLastTrkEtaRelTagVarCSV").c_str()   ,Jet_nLastTrkEtaRelTagVarCSV  );

      tree->SetBranchAddress((name+"TagVarCSV_trackJetPt").c_str()               ,TagVarCSV_trackJetPt              );
      tree->SetBranchAddress((name+"TagVarCSV_jetNTracks").c_str()               ,TagVarCSV_jetNTracks              );
      tree->SetBranchAddress((name+"TagVarCSV_jetNTracksEtaRel").c_str()         ,TagVarCSV_jetNTracksEtaRel        );
      tree->SetBranchAddress((name+"TagVarCSV_trackSumJetEtRatio").c_str()       ,TagVarCSV_trackSumJetEtRatio      );
      tree->SetBranchAddress((name+"TagVarCSV_trackSumJetDeltaR").c_str()        ,TagVarCSV_trackSumJetDeltaR       );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip2dValAboveCharm").c_str()  ,TagVarCSV_trackSip2dValAboveCharm );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip2dSigAboveCharm").c_str()  ,TagVarCSV_trackSip2dSigAboveCharm );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip3dValAboveCharm").c_str()  ,TagVarCSV_trackSip3dValAboveCharm );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip3dSigAboveCharm").c_str()  ,TagVarCSV_trackSip3dSigAboveCharm );
      tree->SetBranchAddress((name+"TagVarCSV_vertexCategory").c_str()           ,TagVarCSV_vertexCategory          );
      tree->SetBranchAddress((name+"TagVarCSV_jetNSecondaryVertices").c_str()    ,TagVarCSV_jetNSecondaryVertices   );
      tree->SetBranchAddress((name+"TagVarCSV_vertexMass").c_str()               ,TagVarCSV_vertexMass              );
      tree->SetBranchAddress((name+"TagVarCSV_vertexNTracks").c_str()            ,TagVarCSV_vertexNTracks           );
      tree->SetBranchAddress((name+"TagVarCSV_vertexEnergyRatio").c_str()        ,TagVarCSV_vertexEnergyRatio       );
      tree->SetBranchAddress((name+"TagVarCSV_vertexJetDeltaR").c_str()          ,TagVarCSV_vertexJetDeltaR         );
      tree->SetBranchAddress((name+"TagVarCSV_flightDistance2dVal").c_str()      ,TagVarCSV_flightDistance2dVal     );
      tree->SetBranchAddress((name+"TagVarCSV_flightDistance2dSig").c_str()      ,TagVarCSV_flightDistance2dSig     );
      tree->SetBranchAddress((name+"TagVarCSV_flightDistance3dVal").c_str()      ,TagVarCSV_flightDistance3dVal     );
      tree->SetBranchAddress((name+"TagVarCSV_flightDistance3dSig").c_str()      ,TagVarCSV_flightDistance3dSig     );

      tree->SetBranchAddress((name+"nTrkTagVarCSV").c_str()               ,&nTrkTagVarCSV             );
      tree->SetBranchAddress((name+"nTrkEtaRelTagVarCSV").c_str()         ,&nTrkEtaRelTagVarCSV       );
      tree->SetBranchAddress((name+"TagVarCSV_trackMomentum").c_str()     ,TagVarCSV_trackMomentum    );
      tree->SetBranchAddress((name+"TagVarCSV_trackEta").c_str()          ,TagVarCSV_trackEta         );
      tree->SetBranchAddress((name+"TagVarCSV_trackPhi").c_str()          ,TagVarCSV_trackPhi         );
      tree->SetBranchAddress((name+"TagVarCSV_trackPtRel").c_str()        ,TagVarCSV_trackPtRel       );
      tree->SetBranchAddress((name+"TagVarCSV_trackPPar").c_str()         ,TagVarCSV_trackPPar        );
      tree->SetBranchAddress((name+"TagVarCSV_trackDeltaR").c_str()       ,TagVarCSV_trackDeltaR      );
      tree->SetBranchAddress((name+"TagVarCSV_trackPtRatio").c_str()      ,TagVarCSV_trackPtRatio     );
      tree->SetBranchAddress((name+"TagVarCSV_trackPParRatio").c_str()    ,TagVarCSV_trackPParRatio   );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip2dVal").c_str()     ,TagVarCSV_trackSip2dVal    );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip2dSig").c_str()     ,TagVarCSV_trackSip2dSig    );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip3dVal").c_str()     ,TagVarCSV_trackSip3dVal    );
      tree->SetBranchAddress((name+"TagVarCSV_trackSip3dSig").c_str()     ,TagVarCSV_trackSip3dSig    );
      tree->SetBranchAddress((name+"TagVarCSV_trackDecayLenVal").c_str()  ,TagVarCSV_trackDecayLenVal );
      tree->SetBranchAddress((name+"TagVarCSV_trackDecayLenSig").c_str()  ,TagVarCSV_trackDecayLenSig );
      tree->SetBranchAddress((name+"TagVarCSV_trackJetDistVal").c_str()   ,TagVarCSV_trackJetDistVal  );
      tree->SetBranchAddress((name+"TagVarCSV_trackJetDistSig").c_str()   ,TagVarCSV_trackJetDistSig  );
      tree->SetBranchAddress((name+"TagVarCSV_trackEtaRel").c_str()       ,TagVarCSV_trackEtaRel      );

    }

    void ReadSubJetSpecificTree(TTree *tree, std::string name="") {
      if(name!="") name += ".";
      tree->SetBranchAddress((name+"Jet_FatJetIdx").c_str(),   Jet_FatJetIdx   );
    }

    void ReadFatJetSpecificTree(TTree *tree, std::string name="") {
      if(name!="") name += ".";
      tree->SetBranchAddress((name+"Jet_ptGroomed").c_str(),   Jet_ptGroomed    );
      tree->SetBranchAddress((name+"Jet_jesGroomed").c_str(),  Jet_jesGroomed   );
      tree->SetBranchAddress((name+"Jet_etaGroomed").c_str(),  Jet_etaGroomed   );
      tree->SetBranchAddress((name+"Jet_phiGroomed").c_str(),  Jet_phiGroomed   );
      tree->SetBranchAddress((name+"Jet_massGroomed").c_str(), Jet_massGroomed  );
      tree->SetBranchAddress((name+"Jet_tau1").c_str(),        Jet_tau1        );
      tree->SetBranchAddress((name+"Jet_tau2").c_str(),        Jet_tau2        );
      tree->SetBranchAddress((name+"Jet_nSubJets").c_str(),    Jet_nSubJets    );
      tree->SetBranchAddress((name+"Jet_nFirstSJ").c_str(),    Jet_nFirstSJ    );
      tree->SetBranchAddress((name+"Jet_nLastSJ").c_str(),     Jet_nLastSJ     );
      tree->SetBranchAddress((name+"nSubJet").c_str(),         &nSubJet        );
      tree->SetBranchAddress((name+"SubJetIdx").c_str(),       SubJetIdx       );
      tree->SetBranchAddress((name+"Jet_nsharedtracks").c_str(), Jet_nsharedtracks );
      tree->SetBranchAddress((name+"Jet_nsubjettracks").c_str(), Jet_nsubjettracks );
      tree->SetBranchAddress((name+"Jet_nsharedsubjettracks").c_str(), Jet_nsharedsubjettracks );
    }
};

#endif
