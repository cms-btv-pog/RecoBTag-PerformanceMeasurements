#ifndef EVENTINFOBRANCHES_H
#define EVENTINFOBRANCHES_H

#include <TTree.h>

const UInt_t nMaxPVs_= 1000;
const UInt_t nMaxPUs_= 1000;
const UInt_t nMaxTrkAll_ = 100000;

class EventInfoBranches {

  public :

    int   nBitTrigger;
    int   BitTrigger[100];
    int   Run;
    int   Evt;
    int   LumiBlock;
    float PVz;
    float PVez;
    float GenPVz;
    float pthat;
    float mcweight;

    int   nPV;
    int   BX;
    float PV_x[nMaxPVs_];
    float PV_y[nMaxPVs_];
    float PV_z[nMaxPVs_];
    float PV_ex[nMaxPVs_];
    float PV_ey[nMaxPVs_];
    float PV_ez[nMaxPVs_];
    float PV_chi2[nMaxPVs_];
    float PV_ndf[nMaxPVs_];
    int   PV_isgood[nMaxPVs_];
    int   PV_isfake[nMaxPVs_];

    float nPUtrue;                 // the true number of pileup interactions that have been added to the event
    int   nPU;                     // the number of pileup interactions that have been added to the event
    int   PU_bunch[nMaxPUs_];      // 0 if on time pileup, -1 or +1 if out-of-time
    float PU_z[nMaxPUs_];          // the true primary vertex position along the z axis for each added interaction
    float PU_sumpT_low[nMaxPUs_];  // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > low_cut
    float PU_sumpT_high[nMaxPUs_]; // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > high_cut
    int   PU_ntrks_low[nMaxPUs_];  // the number of tracks originating from each interaction, where track pT > low_cu
    int   PU_ntrks_high[nMaxPUs_]; // the number of tracks originating from each interaction, where track pT > high_cut

    int   ncQuarks;
    float cQuark_pT[1000];
    float cQuark_eta[1000];
    float cQuark_phi[1000];
    int   cQuark_pdgID[1000];
    int   cQuark_status[1000];
    int   cQuark_fromGSP[1000];

    int   nbQuarks;
    float bQuark_pT[1000];
    float bQuark_eta[1000];
    float bQuark_phi[1000];
    int   bQuark_pdgID[1000];
    int   bQuark_status[1000];
    int   bQuark_fromGSP[1000];

    int   nBHadrons;
    float BHadron_pT[100];
    float BHadron_eta[100];
    float BHadron_phi[100];
    float BHadron_mass[100];
    int   BHadron_pdgID[100];
    int   BHadron_mother[100];
    int   BHadron_hasBdaughter[100];
    float BHadron_SVx[100];
    float BHadron_SVy[100];
    float BHadron_SVz[100];
    int   BHadron_nCharged[100];
    int   BHadron_DHadron1[100];
    int   BHadron_DHadron2[100];

    int   nDHadrons;
    int   nDaughters;
    float DHadron_pT[100];
    float DHadron_eta[100];
    float DHadron_phi[100];
    float DHadron_mass[100];
    float DHadron_SVx[100];
    float DHadron_SVy[100];
    float DHadron_SVz[100];
    int   DHadron_pdgID[100];
    int   DHadron_nDaughters[100];
    int   DHadron_DaughtersPdgID[1500]; // sum(DHadron_nDaughters[i]): needed for daughter pdgIDs
    int   DHadron_nChargedDaughters[100];
    int   DHadron_nCharged[100];

    int   nGenlep;
    float Genlep_pT[100];
    float Genlep_eta[100];
    float Genlep_phi[100];
    int   Genlep_pdgID[100];
    int   Genlep_status[100];
    int   Genlep_mother[100];

    int   nGenquark;
    float Genquark_pT[100];
    float Genquark_eta[100];
    float Genquark_phi[100];
    int   Genquark_pdgID[100];
    int   Genquark_mother[100];

    int   nGenV0;
    float GenV0_pT[1000];
    float GenV0_eta[1000];
    float GenV0_phi[1000];
    int   GenV0_pdgID[1000];
    float GenV0_SVx[1000];
    float GenV0_SVy[1000];
    float GenV0_SVz[1000];
    int   GenV0_nCharged[1000];

    // HP means hard process
    int   nGenPruned;
    float GenPruned_pT[1000];
    float GenPruned_eta[1000];
    float GenPruned_phi[1000];
    float GenPruned_mass[1000];
    int   GenPruned_status[1000];
    int   GenPruned_pdgID[1000];
    int   GenPruned_mother[1000];

    int   nTrkAll;
    float TrkAll_d0[nMaxTrk_];
    float TrkAll_dz[nMaxTrk_];
    float TrkAll_p[nMaxTrk_];
    float TrkAll_pt[nMaxTrk_];
    float TrkAll_eta[nMaxTrk_];
    float TrkAll_phi[nMaxTrk_];
    float TrkAll_chi2[nMaxTrk_];
    int   TrkAll_charge[nMaxTrk_];
    int   TrkAll_nHitAll[nMaxTrk_];
    int   TrkAll_nHitPixel[nMaxTrk_];
    int   TrkAll_nHitStrip[nMaxTrk_];
    int   TrkAll_nHitTIB[nMaxTrk_];
    int   TrkAll_nHitTID[nMaxTrk_];
    int   TrkAll_nHitTOB[nMaxTrk_];
    int   TrkAll_nHitTEC[nMaxTrk_];
    int   TrkAll_nHitPXB[nMaxTrk_];
    int   TrkAll_nHitPXF[nMaxTrk_];
    int   TrkAll_isHitL1[nMaxTrk_];
    int   TrkAll_nSiLayers[nMaxTrk_];
    int   TrkAll_nPxLayers[nMaxTrk_];
    int   TrkAll_PV[nMaxTrk_];
    float TrkAll_PVweight[nMaxTrk_];

    int   nPatMuon;
    int   PatMuon_isGlobal[1000];
    int   PatMuon_isPF[1000];
    int   PatMuon_nTkHit[1000];
    int   PatMuon_nPixHit[1000];
    int   PatMuon_nOutHit[1000];
    int   PatMuon_nMuHit[1000];
    int   PatMuon_nMatched[1000];
    float PatMuon_chi2[1000];
    float PatMuon_chi2Tk[1000];
    float PatMuon_pt[1000];
    float PatMuon_eta[1000];
    float PatMuon_phi[1000];
    float PatMuon_vz[1000];
    float PatMuon_IP[1000];
    float PatMuon_IPsig[1000];
    float PatMuon_IP2D[1000];
    float PatMuon_IP2Dsig[1000];

    int   ttbar_chan, ttbar_trigWord, ttbar_metfilterWord, ttbar_allmepartons, ttbar_matchmepartons;
    int   ttbar_ng,ttbar_gid[25];
    int   ttbar_nl, ttbar_lid[10], ttbar_lgid[10], ttbar_lch[10];
    int   ttbar_nw;
    float ttbar_gpt[25],ttbar_geta[25],ttbar_gphi[25],ttbar_gm[25];
    float ttbar_lpt[10], ttbar_leta[10], ttbar_lphi[10], ttbar_lm[10];
    float ttbar_metpt,ttbar_metphi;
    float ttbar_rho;
    float ttbar_w[500];


    void RegisterTree(TTree *tree) {
      tree->Branch("nBitTrigger", &nBitTrigger,  "nBitTrigger/I");
      tree->Branch("BitTrigger" , BitTrigger  ,  "BitTrigger[nBitTrigger]/I");
      tree->Branch("Run"        , &Run        ,  "Run/I");
      tree->Branch("Evt"        , &Evt        ,  "Evt/I");
      tree->Branch("LumiBlock"  , &LumiBlock  ,  "LumiBlock/I");
      tree->Branch("pthat"      , &pthat      ,  "pthat/F");
      tree->Branch("mcweight"   , &mcweight   ,  "mcweight/F");
      tree->Branch("BX"         , &BX         ,  "BX/I");
      tree->Branch("nPV"        , &nPV        ,  "nPV/I");
      tree->Branch("PVz"        , &PVz        ,  "PVz/F");
      tree->Branch("PVez"       , &PVez       ,  "PVez/F");
      tree->Branch("GenPVz"     , &GenPVz     ,  "GenPVz/F");

      tree->Branch("nPUtrue"      , &nPUtrue     , "nPUtrue/F");
      tree->Branch("nPU"          , &nPU         , "nPU/I"    );
      tree->Branch("PU_bunch"     , PU_bunch     , "PU_bunch[nPU]/I");
      tree->Branch("PU_z"         , PU_z         , "PU_z[nPU]/F");
      tree->Branch("PU_sumpT_low" , PU_sumpT_low , "PU_sumpT_low[nPU]/F");
      tree->Branch("PU_sumpT_high", PU_sumpT_high, "PU_sumpT_high[nPU]/F");
      tree->Branch("PU_ntrks_low" , PU_ntrks_low , "PU_ntrks_low[nPU]/I");
      tree->Branch("PU_ntrks_high", PU_ntrks_high, "PU_ntrks_high[nPU]/I");

      tree->Branch("ncQuarks"    , &ncQuarks       ,"ncQuarks/I");
      tree->Branch("cQuark_pT"     , cQuark_pT     , "cQuark_pT[ncQuarks]/F");
      tree->Branch("cQuark_eta"    , cQuark_eta    , "cQuark_eta[ncQuarks]/F");
      tree->Branch("cQuark_phi"    , cQuark_phi    , "cQuark_phi[ncQuarks]/F");
      tree->Branch("cQuark_pdgID"  , cQuark_pdgID  , "cQuark_pdgID[ncQuarks]/I");
      tree->Branch("cQuark_status" , cQuark_status , "cQuark_status[ncQuarks]/I");
      tree->Branch("cQuark_fromGSP", cQuark_fromGSP, "cQuark_fromGSP[ncQuarks]/I");

      tree->Branch("nbQuarks",        &nbQuarks       ,"nbQuarks/I");
      tree->Branch("bQuark_pT",         bQuark_pT     , "bQuark_pT[nbQuarks]/F");
      tree->Branch("bQuark_eta",        bQuark_eta    , "bQuark_eta[nbQuarks]/F");
      tree->Branch("bQuark_phi",        bQuark_phi    , "bQuark_phi[nbQuarks]/F");
      tree->Branch("bQuark_pdgID",      bQuark_pdgID  , "bQuark_pdgID[nbQuarks]/I");
      tree->Branch("bQuark_status",     bQuark_status , "bQuark_status[nbQuarks]/I");
      tree->Branch("bQuark_fromGSP",    bQuark_fromGSP, "bQuark_fromGSP[nbQuarks]/I");

      tree->Branch("nBHadrons",          &nBHadrons            ,"nBHadrons/I");
      tree->Branch("BHadron_pT",           BHadron_pT          , "BHadron_pT[nBHadrons]/F");
      tree->Branch("BHadron_eta",          BHadron_eta         , "BHadron_eta[nBHadrons]/F");
      tree->Branch("BHadron_phi",          BHadron_phi         , "BHadron_phi[nBHadrons]/F");
      tree->Branch("BHadron_mass",         BHadron_mass        , "BHadron_mass[nBHadrons]/F");
      tree->Branch("BHadron_pdgID",        BHadron_pdgID       , "BHadron_pdgID[nBHadrons]/I");
      tree->Branch("BHadron_mother",       BHadron_mother      , "BHadron_mother[nBHadrons]/I");
      tree->Branch("BHadron_hasBdaughter", BHadron_hasBdaughter, "BHadron_hasBdaughter[nBHadrons]/I");
      tree->Branch("BHadron_SVx",          BHadron_SVx         ,"BHadron_SVx[nBHadrons]/F");
      tree->Branch("BHadron_SVy",          BHadron_SVy         ,"BHadron_SVy[nBHadrons]/F");
      tree->Branch("BHadron_SVz",          BHadron_SVz         ,"BHadron_SVz[nBHadrons]/F");
      tree->Branch("BHadron_nCharged",     BHadron_nCharged    ,"BHadron_nCharged[nBHadrons]/I");
      tree->Branch("BHadron_DHadron1",     BHadron_DHadron1    ,"BHadron_DHadron1[nBHadrons]/I");
      tree->Branch("BHadron_DHadron2",     BHadron_DHadron2    ,"BHadron_DHadron2[nBHadrons]/I");

      tree->Branch("nDHadrons",    &nDHadrons    ,"nDHadrons/I");
      tree->Branch("nDaughters",   &nDaughters   ,"nDaughters/I");
      tree->Branch("DHadron_pT",    DHadron_pT    ,"DHadron_pT[nDHadrons]/F");
      tree->Branch("DHadron_eta",   DHadron_eta   ,"DHadron_eta[nDHadrons]/F");
      tree->Branch("DHadron_phi",   DHadron_phi   ,"DHadron_phi[nDHadrons]/F");
      tree->Branch("DHadron_pdgID", DHadron_pdgID ,"DHadron_pdgID[nDHadrons]/I");
      tree->Branch("DHadron_mass",  DHadron_mass  ,"DHadron_mass[nDHadrons]/F");
      tree->Branch("DHadron_SVx",   DHadron_SVx   ,"DHadron_SVx[nDHadrons]/F");
      tree->Branch("DHadron_SVy",   DHadron_SVy   ,"DHadron_SVy[nDHadrons]/F");
      tree->Branch("DHadron_SVz",   DHadron_SVz   ,"DHadron_SVz[nDHadrons]/F");
      tree->Branch("DHadron_nDaughters",         DHadron_nDaughters        ,"DHadron_nDaughters[nDHadrons]/I");
      tree->Branch("DHadron_DaughtersPdgID",     DHadron_DaughtersPdgID    ,"DHadron_DaughtersPdgID[nDaughters]/I");
      tree->Branch("DHadron_nChargedDaughters",  DHadron_nChargedDaughters ,"DHadron_nChargedDaughters[nDHadrons]/I");
      tree->Branch("DHadron_nCharged", DHadron_nCharged ,"DHadron_nCharged[nDHadrons]/I");

      tree->Branch("nGenlep",     &nGenlep       ,"nGenlep/I");
      tree->Branch("Genlep_pT",     Genlep_pT    , "Genlep_pT[nGenlep]/F");
      tree->Branch("Genlep_eta",    Genlep_eta   , "Genlep_eta[nGenlep]/F");
      tree->Branch("Genlep_phi",    Genlep_phi   , "Genlep_phi[nGenlep]/F");
      tree->Branch("Genlep_pdgID",  Genlep_pdgID , "Genlep_pdgID[nGenlep]/I");
      tree->Branch("Genlep_status", Genlep_status, "Genlep_status[nGenlep]/I");
      tree->Branch("Genlep_mother", Genlep_mother, "Genlep_mother[nGenlep]/I");

      tree->Branch("nGenquark",     &nGenquark       ,"nGenquark/I");
      tree->Branch("Genquark_pT",     Genquark_pT    , "Genquark_pT[nGenquark]/F");
      tree->Branch("Genquark_eta",    Genquark_eta   , "Genquark_eta[nGenquark]/F");
      tree->Branch("Genquark_phi",    Genquark_phi   , "Genquark_phi[nGenquark]/F");
      tree->Branch("Genquark_pdgID",  Genquark_pdgID , "Genquark_pdgID[nGenquark]/I");
      tree->Branch("Genquark_mother", Genquark_mother, "Genquark_mother[nGenquark]/I");

      tree->Branch("nGenPruned",     &nGenPruned       ,"nGenPruned/I");
      tree->Branch("GenPruned_pT",     GenPruned_pT    , "GenPruned_pT[nGenPruned]/F");
      tree->Branch("GenPruned_eta",    GenPruned_eta   , "GenPruned_eta[nGenPruned]/F");
      tree->Branch("GenPruned_phi",    GenPruned_phi   , "GenPruned_phi[nGenPruned]/F");
      tree->Branch("GenPruned_mass",    GenPruned_mass   , "GenPruned_mass[nGenPruned]/F");
      tree->Branch("GenPruned_pdgID",  GenPruned_pdgID , "GenPruned_pdgID[nGenPruned]/I");
      tree->Branch("GenPruned_status", GenPruned_status, "GenPruned_status[nGenPruned]/I");
      tree->Branch("GenPruned_mother", GenPruned_mother, "GenPruned_mother[nGenPruned]/I");

      tree->Branch("nGenV0",        &nGenV0	    ,"nGenV0/I");
      tree->Branch("GenV0_pT",        GenV0_pT       ,"GenV0_pT[nGenV0]/F");
      tree->Branch("GenV0_eta",       GenV0_eta      ,"GenV0_eta[nGenV0]/F");
      tree->Branch("GenV0_phi",       GenV0_phi      ,"GenV0_phi[nGenV0]/F");
      tree->Branch("GenV0_pdgID",     GenV0_pdgID    ,"GenV0_pdgID[nGenV0]/I");
      tree->Branch("GenV0_SVx",       GenV0_SVx      ,"GenV0_SVx[nGenV0]/F");
      tree->Branch("GenV0_SVy",       GenV0_SVy      ,"GenV0_SVy[nGenV0]/F");
      tree->Branch("GenV0_SVz",       GenV0_SVz      ,"GenV0_SVz[nGenV0]/F");
      tree->Branch("GenV0_nCharged",  GenV0_nCharged ,"GenV0_nCharged[nGenV0]/I");
    }

    void RegisterJetTrackTree(TTree *tree) {
      tree->Branch("PV_x"     , PV_x     , "PV_x[nPV]/F");
      tree->Branch("PV_y"     , PV_y     , "PV_y[nPV]/F");
      tree->Branch("PV_z"     , PV_z     , "PV_z[nPV]/F");
      tree->Branch("PV_ex"    , PV_ex    , "PV_ex[nPV]/F");
      tree->Branch("PV_ey"    , PV_ey    , "PV_ey[nPV]/F");
      tree->Branch("PV_ez"    , PV_ez    , "PV_ez[nPV]/F");
      tree->Branch("PV_chi2"  , PV_chi2  , "PV_chi2[nPV]/F");
      tree->Branch("PV_ndf"   , PV_ndf   , "PV_ndf[nPV]/F");
      tree->Branch("PV_isgood", PV_isgood, "PV_isgood[nPV]/I");
      tree->Branch("PV_isfake", PV_isfake, "PV_isfake[nPV]/I");
    }

    void RegisterAllTrackTree(TTree *tree) {
      tree->Branch("nTrkAll",          &nTrkAll,         "nTrkAll/I"                  );
      tree->Branch("TrkAll_d0",        TrkAll_d0,        "TrkAll_d0[nTrkAll]/F"       );
      tree->Branch("TrkAll_dz",        TrkAll_dz,        "TrkAll_dz[nTrkAll]/F"       );
      tree->Branch("TrkAll_p",         TrkAll_p,         "TrkAll_p[nTrkAll]/F"        );
      tree->Branch("TrkAll_pt",        TrkAll_pt,        "TrkAll_pt[nTrkAll]/F"       );
      tree->Branch("TrkAll_eta",       TrkAll_eta,       "TrkAll_eta[nTrkAll]/F"      );
      tree->Branch("TrkAll_phi",       TrkAll_phi,       "TrkAll_phi[nTrkAll]/F"      );
      tree->Branch("TrkAll_chi2",      TrkAll_chi2,      "TrkAll_chi2[nTrkAll]/F"     );
      tree->Branch("TrkAll_charge",    TrkAll_charge,    "TrkAll_charge[nTrkAll]/I"   );
      tree->Branch("TrkAll_nHitAll",   TrkAll_nHitAll,   "TrkAll_nHitAll[nTrkAll]/I"  );
      tree->Branch("TrkAll_nHitPixel", TrkAll_nHitPixel, "TrkAll_nHitPixel[nTrkAll]/I");
      tree->Branch("TrkAll_nHitStrip", TrkAll_nHitStrip, "TrkAll_nHitStrip[nTrkAll]/I");
      tree->Branch("TrkAll_nHitTIB",   TrkAll_nHitTIB,   "TrkAll_nHitTIB[nTrkAll]/I"  );
      tree->Branch("TrkAll_nHitTID",   TrkAll_nHitTID,   "TrkAll_nHitTID[nTrkAll]/I"  );
      tree->Branch("TrkAll_nHitTOB",   TrkAll_nHitTOB,   "TrkAll_nHitTOB[nTrkAll]/I"  );
      tree->Branch("TrkAll_nHitTEC",   TrkAll_nHitTEC,   "TrkAll_nHitTEC[nTrkAll]/I"  );
      tree->Branch("TrkAll_nHitPXB",   TrkAll_nHitPXB,   "TrkAll_nHitPXB[nTrkAll]/I"  );
      tree->Branch("TrkAll_nHitPXF",   TrkAll_nHitPXF,   "TrkAll_nHitPXF[nTrkAll]/I"  );
      tree->Branch("TrkAll_isHitL1",   TrkAll_isHitL1,   "TrkAll_isHitL1[nTrkAll]/I"  );
      tree->Branch("TrkAll_nSiLayers", TrkAll_nSiLayers, "TrkAll_nSiLayers[nTrkAll]/I");
      tree->Branch("TrkAll_nPxLayers", TrkAll_nPxLayers, "TrkAll_nPxLayers[nTrkAll]/I");
      tree->Branch("TrkAll_PV",        TrkAll_PV,        "TrkAll_PV[nTrkAll]/I"       );
      tree->Branch("TrkAll_PVweight",  TrkAll_PVweight,  "TrkAll_PVweight[nTrkAll]/F" );
    }

    void RegisterTTbarTree(TTree *tree) {
      tree->Branch("ttbar_chan" , &ttbar_chan  , "ttbar_chan/I");
      tree->Branch("ttbar_trigWord" , &ttbar_trigWord, "ttbar_trigWord/I");
      tree->Branch("ttbar_metfilterWord" , &ttbar_metfilterWord, "ttbar_metfilterWord/I");
      tree->Branch("ttbar_allmepartons" , &ttbar_allmepartons, "ttbar_allmepartons/I");
      tree->Branch("ttbar_matchmepartons" , &ttbar_matchmepartons, "ttbar_matchmepartons/I");
      tree->Branch("ttbar_ng"   , &ttbar_ng    , "ttbar_ng/I");
      tree->Branch("ttbar_gpt"  ,  ttbar_gpt   , "ttbar_gpt[ttbar_ng]/F");
      tree->Branch("ttbar_geta" ,  ttbar_geta  , "ttbar_geta[ttbar_ng]/F");
      tree->Branch("ttbar_gphi" ,  ttbar_gphi  , "ttbar_gphi[ttbar_ng]/F");
      tree->Branch("ttbar_gm"   ,  ttbar_gm    , "ttbar_gm[ttbar_ng]/F");
      tree->Branch("ttbar_gid"  ,  ttbar_gid   , "ttbar_gid[ttbar_ng]/I");
      tree->Branch("ttbar_nl"   , &ttbar_nl    , "ttbar_nl/I");
      tree->Branch("ttbar_lpt"  ,  ttbar_lpt   , "ttbar_lpt[ttbar_nl]/F");
      tree->Branch("ttbar_leta" ,  ttbar_leta  , "ttbar_leta[ttbar_nl]/F");
      tree->Branch("ttbar_lphi" ,  ttbar_lphi  , "ttbar_lphi[ttbar_nl]/F");
      tree->Branch("ttbar_lm"   ,  ttbar_lm    , "ttbar_lm[ttbar_nl]/F");
      tree->Branch("ttbar_lid"  ,  ttbar_lid   , "ttbar_lid[ttbar_nl]/I");
      tree->Branch("ttbar_lgid" ,  ttbar_lgid  , "ttbar_lgid[ttbar_nl]/I");
      tree->Branch("ttbar_lch"  ,  ttbar_lch   , "ttbar_lch[ttbar_nl]/I");
      tree->Branch("ttbar_metpt",  &ttbar_metpt,  "ttbar_metpt/F");
      tree->Branch("ttbar_metphi", &ttbar_metphi, "ttbar_metphi/F");
      tree->Branch("ttbar_rho"  , &ttbar_rho   , "ttbar_rho/F");
      tree->Branch("ttbar_nw"  ,  &ttbar_nw    , "ttbar_nw/I");
      tree->Branch("ttbar_w"    ,  ttbar_w     , "ttbar_w[ttbar_nw]/F");
    }

    void RegisterPatMuonTree(TTree *tree) {
      tree->Branch("nPatMuon"        , &nPatMuon       , "nPatMuon/I");
      tree->Branch("PatMuon_nMuHit"  , PatMuon_nMuHit  , "PatMuon_nMuHit[nPatMuon]/I");
      tree->Branch("PatMuon_nTkHit"  , PatMuon_nTkHit  , "PatMuon_nTkHit[nPatMuon]/I");
      tree->Branch("PatMuon_nPixHit" , PatMuon_nPixHit , "PatMuon_nPixHit[nPatMuon]/I");
      tree->Branch("PatMuon_nOutHit" , PatMuon_nOutHit , "PatMuon_nOutHit[nPatMuon]/I");
      tree->Branch("PatMuon_isGlobal", PatMuon_isGlobal, "PatMuon_isGlobal[nPatMuon]/I");
      tree->Branch("PatMuon_isPF"    , PatMuon_isPF    , "PatMuon_isPF[nPatMuon]/I");
      tree->Branch("PatMuon_nMatched", PatMuon_nMatched, "PatMuon_nMatched[nPatMuon]/I");
      tree->Branch("PatMuon_chi2"    , PatMuon_chi2    , "PatMuon_chi2[nPatMuon]/F");
      tree->Branch("PatMuon_chi2Tk"  , PatMuon_chi2Tk  , "PatMuon_chi2Tk[nPatMuon]/F");
      tree->Branch("PatMuon_pt"      , PatMuon_pt      , "PatMuon_pt[nPatMuon]/F");
      tree->Branch("PatMuon_eta"     , PatMuon_eta     , "PatMuon_eta[nPatMuon]/F");
      tree->Branch("PatMuon_phi"     , PatMuon_phi     , "PatMuon_phi[nPatMuon]/F");
      tree->Branch("PatMuon_vz"      , PatMuon_vz      , "PatMuon_vz[nPatMuon]/F");
      tree->Branch("PatMuon_IP"      , PatMuon_IP      , "PatMuon_IP[nPatMuon]/F");
      tree->Branch("PatMuon_IPsig"   , PatMuon_IPsig   , "PatMuon_IPsig[nPatMuon]/F");
      tree->Branch("PatMuon_IP2D"    , PatMuon_IP2D    , "PatMuon_IP2D[nPatMuon]/F");
      tree->Branch("PatMuon_IP2Dsig" , PatMuon_IP2Dsig , "PatMuon_IP2Dsig[nPatMuon]/F");
    }

    //------------------------------------------------------------------------------------------------------------------

    void ReadTree(TTree *tree) {
      tree->SetBranchAddress("nBitTrigger", &nBitTrigger);
      tree->SetBranchAddress("BitTrigger" , BitTrigger  );
      tree->SetBranchAddress("Run"        , &Run        );
      tree->SetBranchAddress("Evt"        , &Evt        );
      tree->SetBranchAddress("LumiBlock"  , &LumiBlock  );
      tree->SetBranchAddress("pthat"      , &pthat      );
      tree->SetBranchAddress("mcweight"   , &mcweight   );
      tree->SetBranchAddress("nPV"        , &nPV        );
      tree->SetBranchAddress("PVz"        , &PVz        );
      tree->SetBranchAddress("PVez"       , &PVez       );
      tree->SetBranchAddress("GenPVz"     , &GenPVz     );

      tree->SetBranchAddress("nPUtrue"      , &nPUtrue     );
      tree->SetBranchAddress("nPU"          , &nPU         );
      tree->SetBranchAddress("PU_bunch"     , PU_bunch     );
      tree->SetBranchAddress("PU_z"         , PU_z         );
      tree->SetBranchAddress("PU_sumpT_low" , PU_sumpT_low );
      tree->SetBranchAddress("PU_sumpT_high", PU_sumpT_high);
      tree->SetBranchAddress("PU_ntrks_low" , PU_ntrks_low );
      tree->SetBranchAddress("PU_ntrks_high", PU_ntrks_high);

      tree->SetBranchAddress("ncQuarks"      , &ncQuarks     );
      tree->SetBranchAddress("cQuark_pT"     , cQuark_pT     );
      tree->SetBranchAddress("cQuark_eta"    , cQuark_eta    );
      tree->SetBranchAddress("cQuark_phi"    , cQuark_phi    );
      tree->SetBranchAddress("cQuark_pdgID"  , cQuark_pdgID  );
      tree->SetBranchAddress("cQuark_status" , cQuark_status );
      tree->SetBranchAddress("cQuark_fromGSP", cQuark_fromGSP);

      tree->SetBranchAddress("nbQuarks",          &nbQuarks     );
      tree->SetBranchAddress("bQuark_pT",         bQuark_pT     );
      tree->SetBranchAddress("bQuark_eta",        bQuark_eta    );
      tree->SetBranchAddress("bQuark_phi",        bQuark_phi    );
      tree->SetBranchAddress("bQuark_pdgID",      bQuark_pdgID  );
      tree->SetBranchAddress("bQuark_status",     bQuark_status );
      tree->SetBranchAddress("bQuark_fromGSP",    bQuark_fromGSP);

      tree->SetBranchAddress("nBHadrons",            &nBHadrons          );
      tree->SetBranchAddress("BHadron_pT",           BHadron_pT          );
      tree->SetBranchAddress("BHadron_eta",          BHadron_eta         );
      tree->SetBranchAddress("BHadron_phi",          BHadron_phi         );
      tree->SetBranchAddress("BHadron_mass",         BHadron_mass        );
      tree->SetBranchAddress("BHadron_pdgID",        BHadron_pdgID       );
      tree->SetBranchAddress("BHadron_mother",       BHadron_mother      );
      tree->SetBranchAddress("BHadron_hasBdaughter", BHadron_hasBdaughter);
      tree->SetBranchAddress("BHadron_SVx",          BHadron_SVx);
      tree->SetBranchAddress("BHadron_SVy",          BHadron_SVy);
      tree->SetBranchAddress("BHadron_SVz",          BHadron_SVz);
      tree->SetBranchAddress("BHadron_nCharged",     BHadron_nCharged);
      tree->SetBranchAddress("BHadron_DHadron1",     BHadron_DHadron1);
      tree->SetBranchAddress("BHadron_DHadron2",     BHadron_DHadron2);

      tree->SetBranchAddress("nDHadrons",    &nDHadrons    );
      tree->SetBranchAddress("nDaughters",   &nDaughters   );
      tree->SetBranchAddress("DHadron_pT",    DHadron_pT    );
      tree->SetBranchAddress("DHadron_eta",   DHadron_eta   );
      tree->SetBranchAddress("DHadron_phi",   DHadron_phi   );
      tree->SetBranchAddress("DHadron_pdgID", DHadron_pdgID );
      tree->SetBranchAddress("DHadron_mass",  DHadron_mass  );
      tree->SetBranchAddress("DHadron_SVx",   DHadron_SVx   );
      tree->SetBranchAddress("DHadron_SVy",   DHadron_SVy   );
      tree->SetBranchAddress("DHadron_SVz",   DHadron_SVz   );
      tree->SetBranchAddress("DHadron_nDaughters",         DHadron_nDaughters        );
      tree->SetBranchAddress("DHadron_DaughtersPdgID",     DHadron_DaughtersPdgID    );
      tree->SetBranchAddress("DHadron_nChargedDaughters",  DHadron_nChargedDaughters );
      tree->SetBranchAddress("DHadron_nCharged", DHadron_nCharged );

      tree->SetBranchAddress("nGenlep",       &nGenlep     );
      tree->SetBranchAddress("Genlep_pT",     Genlep_pT    );
      tree->SetBranchAddress("Genlep_eta",    Genlep_eta   );
      tree->SetBranchAddress("Genlep_phi",    Genlep_phi   );
      tree->SetBranchAddress("Genlep_pdgID",  Genlep_pdgID );
      tree->SetBranchAddress("Genlep_status", Genlep_status);
      tree->SetBranchAddress("Genlep_mother", Genlep_mother);

      tree->SetBranchAddress("nGenquark",       &nGenquark     );
      tree->SetBranchAddress("Genquark_pT",     Genquark_pT    );
      tree->SetBranchAddress("Genquark_eta",    Genquark_eta   );
      tree->SetBranchAddress("Genquark_phi",    Genquark_phi   );
      tree->SetBranchAddress("Genquark_pdgID",  Genquark_pdgID );
      tree->SetBranchAddress("Genquark_mother", Genquark_mother);

      tree->SetBranchAddress("nGenPruned",       &nGenPruned     );
      tree->SetBranchAddress("GenPruned_pT",     GenPruned_pT    );
      tree->SetBranchAddress("GenPruned_eta",    GenPruned_eta   );
      tree->SetBranchAddress("GenPruned_phi",    GenPruned_phi   );
      tree->SetBranchAddress("GenPruned_mass",    GenPruned_mass   );
      tree->SetBranchAddress("GenPruned_pdgID",  GenPruned_pdgID );
      tree->SetBranchAddress("GenPruned_status", GenPruned_status);
      tree->SetBranchAddress("GenPruned_mother", GenPruned_mother);

      tree->SetBranchAddress("nGenV0",       &nGenV0 );
      tree->SetBranchAddress("GenV0_pT",       GenV0_pT );
      tree->SetBranchAddress("GenV0_eta",      GenV0_eta );
      tree->SetBranchAddress("GenV0_phi",      GenV0_phi );
      tree->SetBranchAddress("GenV0_pdgID",    GenV0_pdgID );
      tree->SetBranchAddress("GenV0_SVx",      GenV0_SVx );
      tree->SetBranchAddress("GenV0_SVy",      GenV0_SVy );
      tree->SetBranchAddress("GenV0_SVz",      GenV0_SVz );
      tree->SetBranchAddress("GenV0_nCharged", GenV0_nCharged );
    }

    void ReadJetTrackTree(TTree *tree) {
      tree->SetBranchAddress("PV_x"     , PV_x     );
      tree->SetBranchAddress("PV_y"     , PV_y     );
      tree->SetBranchAddress("PV_z"     , PV_z     );
      tree->SetBranchAddress("PV_ex"    , PV_ex    );
      tree->SetBranchAddress("PV_ey"    , PV_ey    );
      tree->SetBranchAddress("PV_ez"    , PV_ez    );
      tree->SetBranchAddress("PV_chi2"  , PV_chi2  );
      tree->SetBranchAddress("PV_ndf"   , PV_ndf   );
      tree->SetBranchAddress("PV_isgood", PV_isgood);
      tree->SetBranchAddress("PV_isfake", PV_isfake);
    }

   void ReadAllTrackTree(TTree *tree) {
      tree->SetBranchAddress("nTrkAll",          &nTrkAll         );
      tree->SetBranchAddress("TrkAll_d0",        TrkAll_d0        );
      tree->SetBranchAddress("TrkAll_dz",        TrkAll_dz        );
      tree->SetBranchAddress("TrkAll_p",         TrkAll_p         );
      tree->SetBranchAddress("TrkAll_pt",        TrkAll_pt        );
      tree->SetBranchAddress("TrkAll_eta",       TrkAll_eta       );
      tree->SetBranchAddress("TrkAll_phi",       TrkAll_phi       );
      tree->SetBranchAddress("TrkAll_chi2",      TrkAll_chi2      );
      tree->SetBranchAddress("TrkAll_charge",    TrkAll_charge    );
      tree->SetBranchAddress("TrkAll_nHitAll",   TrkAll_nHitAll   );
      tree->SetBranchAddress("TrkAll_nHitPixel", TrkAll_nHitPixel );
      tree->SetBranchAddress("TrkAll_nHitStrip", TrkAll_nHitStrip );
      tree->SetBranchAddress("TrkAll_nHitTIB",   TrkAll_nHitTIB   );
      tree->SetBranchAddress("TrkAll_nHitTID",   TrkAll_nHitTID   );
      tree->SetBranchAddress("TrkAll_nHitTOB",   TrkAll_nHitTOB   );
      tree->SetBranchAddress("TrkAll_nHitTEC",   TrkAll_nHitTEC   );
      tree->SetBranchAddress("TrkAll_nHitPXB",   TrkAll_nHitPXB   );
      tree->SetBranchAddress("TrkAll_nHitPXF",   TrkAll_nHitPXF   );
      tree->SetBranchAddress("TrkAll_isHitL1",   TrkAll_isHitL1   );
      tree->SetBranchAddress("TrkAll_nSiLayers", TrkAll_nSiLayers );
      tree->SetBranchAddress("TrkAll_nPxLayers", TrkAll_nPxLayers );
      tree->SetBranchAddress("TrkAll_PV",        TrkAll_PV        );
      tree->SetBranchAddress("TrkAll_PVweight",  TrkAll_PVweight  );
    }

    void ReadTTbarTree(TTree *tree) {
      tree->SetBranchAddress("ttbar_chan" , &ttbar_chan);
      tree->SetBranchAddress("ttbar_trigWord" , &ttbar_trigWord);
      tree->SetBranchAddress("ttbar_ng"   , &ttbar_ng);
      tree->SetBranchAddress("ttbar_gpt"  ,  ttbar_gpt);
      tree->SetBranchAddress("ttbar_geta" ,  ttbar_geta);
      tree->SetBranchAddress("ttbar_gphi" ,  ttbar_gphi);
      tree->SetBranchAddress("ttbar_gm"   ,  ttbar_gm);
      tree->SetBranchAddress("ttbar_gid"  ,  ttbar_gid);
      tree->SetBranchAddress("ttbar_nl"   , &ttbar_nl);
      tree->SetBranchAddress("ttbar_lpt"  ,  ttbar_lpt); 
      tree->SetBranchAddress("ttbar_leta" ,  ttbar_leta);
      tree->SetBranchAddress("ttbar_lphi" ,  ttbar_lphi);
      tree->SetBranchAddress("ttbar_lm"   ,  ttbar_lm);
      tree->SetBranchAddress("ttbar_lid"  ,  ttbar_lid);
      tree->SetBranchAddress("ttbar_lgid" ,  ttbar_lgid);
      tree->SetBranchAddress("ttbar_lch"  ,  ttbar_lch);
      tree->SetBranchAddress("ttbar_metpt"  , &ttbar_metpt);
      tree->SetBranchAddress("ttbar_metphi"  , &ttbar_metphi);
      tree->SetBranchAddress("ttbar_rho"  , &ttbar_rho);
      tree->SetBranchAddress("ttbar_nw"  ,  &ttbar_nw);
      tree->SetBranchAddress("ttbar_w"    ,  ttbar_w);
    }

    void ReadPatMuonTree(TTree *tree) {
      tree->SetBranchAddress("nPatMuon"        , &nPatMuon       );
      tree->SetBranchAddress("PatMuon_nMuHit"  , PatMuon_nMuHit  );
      tree->SetBranchAddress("PatMuon_nTkHit"  , PatMuon_nTkHit  );
      tree->SetBranchAddress("PatMuon_nPixHit" , PatMuon_nPixHit );
      tree->SetBranchAddress("PatMuon_nOutHit" , PatMuon_nOutHit );
      tree->SetBranchAddress("PatMuon_isGlobal", PatMuon_isGlobal);
      tree->SetBranchAddress("PatMuon_isPF"    , PatMuon_isPF    );
      tree->SetBranchAddress("PatMuon_nMatched", PatMuon_nMatched);
      tree->SetBranchAddress("PatMuon_chi2"    , PatMuon_chi2    );
      tree->SetBranchAddress("PatMuon_chi2Tk"  , PatMuon_chi2Tk  );
      tree->SetBranchAddress("PatMuon_pt"      , PatMuon_pt      );
      tree->SetBranchAddress("PatMuon_eta"     , PatMuon_eta     );
      tree->SetBranchAddress("PatMuon_phi"     , PatMuon_phi     );
      tree->SetBranchAddress("PatMuon_vz"      , PatMuon_vz      );
      tree->SetBranchAddress("PatMuon_IP"      , PatMuon_IP      );
      tree->SetBranchAddress("PatMuon_IPsig"   , PatMuon_IPsig   );
      tree->SetBranchAddress("PatMuon_IP2D"    , PatMuon_IP2D    );
      tree->SetBranchAddress("PatMuon_IP2Dsig" , PatMuon_IP2Dsig );
    }
};

#endif

