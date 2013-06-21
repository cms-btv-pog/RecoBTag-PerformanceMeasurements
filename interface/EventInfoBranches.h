#ifndef EVENTINFOBRANCHES_H
#define EVENTINFOBRANCHES_H

#include <TTree.h>

const UInt_t nMaxPVs_= 1000;
const UInt_t nMaxPUs_= 1000;

class EventInfoBranches {

  public :

    int   nBitTrigger;
    int   BitTrigger[100];
    int   Run;
    int   Evt;
    int   LumiBlock;
    float PVz;
    float pthat;
    float mcweight;

    int   nPV;
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
    float BHadron_pT[1000];
    float BHadron_eta[1000];
    float BHadron_phi[1000];
    float BHadron_mass[1000];
    int   BHadron_pdgID[1000];
    int   BHadron_status[1000];
    int   BHadron_mother[1000];
    int   BHadron_hasBdaughter[1000];
    
//new
    int   nDHadrons;
    int   nDaughters;    
    float DHadron_pT[100];
    float DHadron_eta[100];
    float DHadron_phi[100];
    float DHadron_mass[100];
    float DHadron_vx[100];
    float DHadron_vy[100];
    float DHadron_vz[100];
    float DHadron_daughterVx[100];
    float DHadron_daughterVy[100];
    float DHadron_daughterVz[100];
    int   DHadron_pdgID[100];
    int   DHadron_nDaughters[100];
    //sum(DHadron_nDaughters[i]): needed for daughter pdgIDs
    int   DHadron_DaughtersPdgID[1500];
    int   DHadron_nChargedDaughters[100];
//new    

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

    int   ttbar_chan;
    float lepton1_pT;
    float lepton2_pT;
    float lepton1_eta;
    float lepton2_eta;
    float lepton1_phi;
    float lepton2_phi;
    float met;
    float mll;
    int   trig_ttbar;


    void RegisterTree(TTree *tree) {
      tree->Branch("nBitTrigger", &nBitTrigger,  "nBitTrigger/I");
      tree->Branch("BitTrigger" , BitTrigger  ,  "BitTrigger[nBitTrigger]/I");
      tree->Branch("Run"        , &Run        ,  "Run/I");
      tree->Branch("Evt"        , &Evt        ,  "Evt/I");
      tree->Branch("LumiBlock"  , &LumiBlock  ,  "LumiBlock/I");
      tree->Branch("pthat"      , &pthat      ,  "pthat/F");
      tree->Branch("mcweight"   , &mcweight   ,  "mcweight/F");
      tree->Branch("nPV"        , &nPV        ,  "nPV/I");
      tree->Branch("PVz"        , &PVz        ,  "PVz/F");

      tree->Branch("nPUtrue"      , &nPUtrue     , "nPUtrue/F");
      tree->Branch("nPU"          , &nPU         , "nPU/I"    );
      // tree->Branch("PU_bunch"     , PU_bunch     , "PU_bunch[nPU]/I");
      // tree->Branch("PU_z"         , PU_z         , "PU_z[nPU]/F");
      // tree->Branch("PU_sumpT_low" , PU_sumpT_low , "PU_sumpT_low[nPU]/F");
      // tree->Branch("PU_sumpT_high", PU_sumpT_high, "PU_sumpT_high[nPU]/F");
      // tree->Branch("PU_ntrks_low" , PU_ntrks_low , "PU_ntrks_low[nPU]/I");
      // tree->Branch("PU_ntrks_high", PU_ntrks_high, "PU_ntrks_high[nPU]/I");

      tree->Branch("ncQuarks"      , &ncQuarks     , "ncQuarks/I");
      tree->Branch("cQuark_pT"     , cQuark_pT     , "cQuark_pT[ncQuarks]/F");
      tree->Branch("cQuark_eta"    , cQuark_eta    , "cQuark_eta[ncQuarks]/F");
      tree->Branch("cQuark_phi"    , cQuark_phi    , "cQuark_phi[ncQuarks]/F");
      tree->Branch("cQuark_pdgID"  , cQuark_pdgID  , "cQuark_pdgID[ncQuarks]/I");
      tree->Branch("cQuark_status" , cQuark_status , "cQuark_status[ncQuarks]/I");
      tree->Branch("cQuark_fromGSP", cQuark_fromGSP, "cQuark_fromGSP[ncQuarks]/I");

      tree->Branch("nbQuarks",          &nbQuarks     , "nbQuarks/I");
      tree->Branch("bQuark_pT",         bQuark_pT     , "bQuark_pT[nbQuarks]/F");
      tree->Branch("bQuark_eta",        bQuark_eta    , "bQuark_eta[nbQuarks]/F");
      tree->Branch("bQuark_phi",        bQuark_phi    , "bQuark_phi[nbQuarks]/F");
      tree->Branch("bQuark_pdgID",      bQuark_pdgID  , "bQuark_pdgID[nbQuarks]/I");
      tree->Branch("bQuark_status",     bQuark_status , "bQuark_status[nbQuarks]/I");
      tree->Branch("bQuark_fromGSP",    bQuark_fromGSP, "bQuark_fromGSP[nbQuarks]/I");

      tree->Branch("nBHadrons",            &nBHadrons          , "nBHadrons/I");
      tree->Branch("BHadron_pT",           BHadron_pT          , "BHadron_pT[nBHadrons]/F");
      tree->Branch("BHadron_eta",          BHadron_eta         , "BHadron_eta[nBHadrons]/F");
      tree->Branch("BHadron_phi",          BHadron_phi         , "BHadron_phi[nBHadrons]/F");
      tree->Branch("BHadron_mass",         BHadron_mass        , "BHadron_mass[nBHadrons]/F");
      tree->Branch("BHadron_pdgID",        BHadron_pdgID       , "BHadron_pdgID[nBHadrons]/I");
      tree->Branch("BHadron_status",       BHadron_status      , "BHadron_status[nBHadrons]/I");
      tree->Branch("BHadron_mother",       BHadron_mother      , "BHadron_mother[nBHadrons]/I");
      tree->Branch("BHadron_hasBdaughter", BHadron_hasBdaughter, "BHadron_hasBdaughter[nBHadrons]/I");
      
//new
      tree->Branch("nDHadrons",            &nDHadrons          , "nDHadrons/I");
      tree->Branch("nDaughters",    &nDaughters   ,"nDaughters/I");
      tree->Branch("DHadron_pT",    DHadron_pT    ,"DHadron_pT[nDHadrons]/F");
      tree->Branch("DHadron_eta",   DHadron_eta   ,"DHadron_eta[nDHadrons]/F");
      tree->Branch("DHadron_phi",   DHadron_phi   ,"DHadron_phi[nDHadrons]/F");
      tree->Branch("DHadron_mass",  DHadron_mass  ,"DHadron_mass[nDHadrons]/F");
      tree->Branch("DHadron_vx",    DHadron_vx    ,"DHadron_vx[nDHadrons]/F");
      tree->Branch("DHadron_vy",    DHadron_vy    ,"DHadron_vy[nDHadrons]/F");
      tree->Branch("DHadron_vz",    DHadron_vz    ,"DHadron_vz[nDHadrons]/F");
      tree->Branch("DHadron_daughterVx",  DHadron_daughterVx  ,"DHadron_daughterVx[nDHadrons]/F");
      tree->Branch("DHadron_daughterVy",  DHadron_daughterVy  ,"DHadron_daughterVy[nDHadrons]/F");
      tree->Branch("DHadron_daughterVz",  DHadron_daughterVz  ,"DHadron_daughterVz[nDHadrons]/F");
      tree->Branch("DHadron_pdgID", DHadron_pdgID ,"DHadron_pdgID[nDHadrons]/I");
      tree->Branch("DHadron_nDaughters", DHadron_nDaughters ,"DHadron_nDaughters[nDHadrons]/I");
      tree->Branch("DHadron_DaughtersPdgID", DHadron_DaughtersPdgID,"DHadron_DaughtersPdgID[nDaughters]/I");
      tree->Branch("DHadron_nChargedDaughters",     &DHadron_nChargedDaughters    ,"DHadron_nChargedDaughters[nDHadrons]/I");
//new      

      tree->Branch("nGenlep",       &nGenlep     , "nGenlep/I");
      tree->Branch("Genlep_pT",     Genlep_pT    , "Genlep_pT[nGenlep]/F");
      tree->Branch("Genlep_eta",    Genlep_eta   , "Genlep_eta[nGenlep]/F");
      tree->Branch("Genlep_phi",    Genlep_phi   , "Genlep_phi[nGenlep]/F");
      tree->Branch("Genlep_pdgID",  Genlep_pdgID , "Genlep_pdgID[nGenlep]/I");
      tree->Branch("Genlep_status", Genlep_status, "Genlep_status[nGenlep]/I");
      tree->Branch("Genlep_mother", Genlep_mother, "Genlep_mother[nGenlep]/I");

      tree->Branch("nGenquark",       &nGenquark     , "nGenquark/I");
      tree->Branch("Genquark_pT",     Genquark_pT    , "Genquark_pT[nGenquark]/F");
      tree->Branch("Genquark_eta",    Genquark_eta   , "Genquark_eta[nGenquark]/F");
      tree->Branch("Genquark_phi",    Genquark_phi   , "Genquark_phi[nGenquark]/F");
      tree->Branch("Genquark_pdgID",  Genquark_pdgID , "Genquark_pdgID[nGenquark]/I");
      tree->Branch("Genquark_mother", Genquark_mother, "Genquark_mother[nGenquark]/I");

    }

    void RegisterJPTree(TTree *tree) {
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

    void RegisterTTbarTree(TTree *tree) {
      tree->Branch("ttbar_chan" , &ttbar_chan , "ttbar_chan/I");
      tree->Branch("lepton1_pT" , &lepton1_pT , "lepton1_pT/F");
      tree->Branch("lepton1_eta", &lepton1_eta, "lepton1_eta/F");
      tree->Branch("lepton1_phi", &lepton1_phi, "lepton1_phi/F");
      tree->Branch("lepton2_pT" , &lepton2_pT , "lepton2_pT/F");
      tree->Branch("lepton2_eta", &lepton2_eta, "lepton2_eta/F");
      tree->Branch("lepton2_phi", &lepton2_phi, "lepton2_phi/F");
      tree->Branch("met"        , &met        , "met/F");
      tree->Branch("mll"        , &mll        , "mll/F");
      tree->Branch("trig_ttbar" , &trig_ttbar , "trig_ttbar/I");
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

      tree->SetBranchAddress("nPUtrue"      , &nPUtrue     );
      tree->SetBranchAddress("nPU"          , &nPU         );
      // tree->SetBranchAddress("PU_bunch"     , PU_bunch     );
      // tree->SetBranchAddress("PU_z"         , PU_z         );
      // tree->SetBranchAddress("PU_sumpT_low" , PU_sumpT_low );
      // tree->SetBranchAddress("PU_sumpT_high", PU_sumpT_high);
      // tree->SetBranchAddress("PU_ntrks_low" , PU_ntrks_low );
      // tree->SetBranchAddress("PU_ntrks_high", PU_ntrks_high);

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
      tree->SetBranchAddress("BHadron_status",       BHadron_status      );
      tree->SetBranchAddress("BHadron_mother",       BHadron_mother      );
      tree->SetBranchAddress("BHadron_hasBdaughter", BHadron_hasBdaughter);

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
    }

    void ReadJPTree(TTree *tree) {
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

   void ReadTTbarTree(TTree *tree) {
      tree->SetBranchAddress("ttbar_chan" , &ttbar_chan );
      tree->SetBranchAddress("lepton1_pT" , &lepton1_pT );
      tree->SetBranchAddress("lepton1_eta", &lepton1_eta);
      tree->SetBranchAddress("lepton1_phi", &lepton1_phi);
      tree->SetBranchAddress("lepton2_pT" , &lepton2_pT );
      tree->SetBranchAddress("lepton2_eta", &lepton2_eta);
      tree->SetBranchAddress("lepton2_phi", &lepton2_phi);
      tree->SetBranchAddress("met"        , &met        );
      tree->SetBranchAddress("mll"        , &mll        );
      tree->SetBranchAddress("trig_ttbar" , &trig_ttbar );
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

