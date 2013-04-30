#ifndef JETINFOBRANCHES_H
#define JETINFOBRANCHES_H

#include <TTree.h>

const UInt_t nMaxJets_ = 10000;

class JetInfoBranches { 

	public : 

  int nJet;
  float Jet_pt[nMaxJets_];
  float Jet_genpt[nMaxJets_];
  float Jet_residual[nMaxJets_];
  float Jet_jes[nMaxJets_];
  float Jet_eta[nMaxJets_];
  float Jet_phi[nMaxJets_];
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
  float Jet_Bprob[nMaxJets_];
  float Jet_BprobP[nMaxJets_];
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

  void RegisterTree(TTree *tree, std::string name="JetInfo") {
    tree->Branch((name+".nJet").c_str(),            &nJet           ,(name+".nJet/I").c_str());
    tree->Branch((name+".Jet_pt").c_str(),          Jet_pt	    ,(name+".Jet_pt["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_genpt").c_str(),       Jet_genpt       ,(name+".Jet_genpt["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_residual").c_str(),    Jet_residual    ,(name+".Jet_residual["+name+".nJet]/F").c_str()); 
    tree->Branch((name+".Jet_jes").c_str(),         Jet_jes         ,(name+".Jet_jes["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_eta").c_str(),         Jet_eta         ,(name+".Jet_eta["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_phi").c_str(),         Jet_phi         ,(name+".Jet_phi["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_ntracks").c_str(),     Jet_ntracks     ,(name+".Jet_ntracks["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_flavour").c_str(),     Jet_flavour     ,(name+".Jet_flavour["+name+".nJet]/I").c_str());
    //  tree->Branch(name+".Jet_Ip1N",	   Jet_Ip1N	   ,(name+".Jet_Ip1N["+name+".nJet]/F").c_str());
    //  tree->Branch(name+".Jet_Ip1P",	   Jet_Ip1P	   ,(name+".Jet_Ip1P["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_Ip2N").c_str(),        Jet_Ip2N        ,(name+".Jet_Ip2N["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_Ip2P").c_str(),        Jet_Ip2P        ,(name+".Jet_Ip2P["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_Ip3N").c_str(),        Jet_Ip3N        ,(name+".Jet_Ip3N["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_Ip3P").c_str(),        Jet_Ip3P        ,(name+".Jet_Ip3P["+name+".nJet]/F").c_str());
    //  tree->Branch(name+".Jet_Ip4N",	   Jet_Ip4N	   ,(name+".Jet_Ip4N["+name+".nJet]/F").c_str());
    //  tree->Branch(name+".Jet_Ip4P",	   Jet_Ip4P	   ,(name+".Jet_Ip4P["+name+".nJet]/F").c_str());
    //  tree->Branch(name+".Jet_Mass4N",	   Jet_Mass4N	   ,(name+".Jet_Mass4N["+name+".nJet]/F").c_str());
    //  tree->Branch(name+".Jet_Mass4P",	   Jet_Mass4P	   ,(name+".Jet_Mass4P["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_ProbaN").c_str(),      Jet_ProbaN     ,(name+".Jet_ProbaN["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_ProbaP").c_str(),      Jet_ProbaP     ,(name+".Jet_ProbaP["+name+".nJet]/F").c_str());
    //$$  tree->Branch((name+".Jet_Proba").c_str(),       Jet_Proba       ,(name+".Jet_Proba["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_BprobN").c_str(),      Jet_BprobN     ,(name+".Jet_BprobN["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_BprobP").c_str(),      Jet_BprobP     ,(name+".Jet_BprobP["+name+".nJet]/F").c_str());
    //$$  tree->Branch((name+".Jet_Bprob").c_str(),       Jet_Bprob       ,(name+".Jet_Bprob["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SvxN").c_str(),        Jet_SvxN       ,(name+".Jet_SvxN["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_Svx").c_str(),         Jet_Svx        ,(name+".Jet_Svx["+name+".nJet]/F").c_str());
    //  tree->Branch((name+".Jet_SvxNTracks").c_str(),  Jet_SvxNTracks  ,(name+".Jet_SvxNTracks["+name+".nJet]/I").c_str());
    //  tree->Branch((name+".Jet_SvxTracks").c_str(),   Jet_SvxTracks   ,(name+".Jet_SvxTracks["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_SvxNHP").c_str(),      Jet_SvxNHP     ,(name+".Jet_SvxNHP["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SvxHP").c_str(),       Jet_SvxHP      ,(name+".Jet_SvxHP["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SvxMass").c_str(),     Jet_SvxMass    ,(name+".Jet_SvxMass["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombSvxN").c_str(),    Jet_CombSvxN   ,(name+".Jet_CombSvxN["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombSvxP").c_str(),    Jet_CombSvxP   ,(name+".Jet_CombSvxP["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombSvx").c_str(),     Jet_CombSvx    ,(name+".Jet_CombSvx["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_RetCombSvxN").c_str(), Jet_RetCombSvxN   ,(name+".Jet_RetCombSvxN["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_RetCombSvxP").c_str(), Jet_RetCombSvxP   ,(name+".Jet_RetCombSvxP["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_RetCombSvx").c_str(), Jet_RetCombSvx    ,(name+".Jet_RetCombSvx["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_CombCSVJP_N").c_str(), Jet_CombCSVJP_N   ,(name+".Jet_CombCSVJP_N["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombCSVJP_P").c_str(), Jet_CombCSVJP_P   ,(name+".Jet_CombCSVJP_P["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombCSVJP").c_str(),Jet_CombCSVJP      ,(name+".Jet_CombCSVJP["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_CombCSVSL_N").c_str(), Jet_CombCSVSL_N    ,(name+".Jet_CombCSVSL_N["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombCSVSL_P").c_str(), Jet_CombCSVSL_P    ,(name+".Jet_CombCSVSL_P["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombCSVSL").c_str(),   Jet_CombCSVSL ,(name+".Jet_CombCSVSL["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_CombCSVJPSL_N").c_str(), Jet_CombCSVJPSL_N  ,(name+".Jet_CombCSVJPSL_N["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombCSVJPSL_P").c_str(), Jet_CombCSVJPSL_P  ,(name+".Jet_CombCSVJPSL_P["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombCSVJPSL").c_str(), Jet_CombCSVJPSL    ,(name+".Jet_CombCSVJPSL["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_SimpIVF_HP").c_str(),  Jet_SimpIVF_HP  ,(name+".Jet_SimpIVF_HP["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SimpIVF_HE").c_str(),  Jet_SimpIVF_HE  ,(name+".Jet_SimpIVF_HE["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_DoubIVF_HE").c_str(),  Jet_DoubIVF_HE  ,(name+".Jet_DoubIVF_HE["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombIVF").c_str(),     Jet_CombIVF     ,(name+".Jet_CombIVF["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_CombIVF_P").c_str(), Jet_CombIVF_P   ,(name+".Jet_CombIVF_P["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_SoftMuN").c_str(),     Jet_SoftMuN     ,(name+".Jet_SoftMuN["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SoftMuP").c_str(),     Jet_SoftMuP     ,(name+".Jet_SoftMuP["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SoftMu").c_str(),      Jet_SoftMu      ,(name+".Jet_SoftMu["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_SoftElN").c_str(),     Jet_SoftElN     ,(name+".Jet_SoftElN["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SoftElP").c_str(),     Jet_SoftElP     ,(name+".Jet_SoftElP["+name+".nJet]/F").c_str());
    tree->Branch((name+".Jet_SoftEl").c_str(),      Jet_SoftEl      ,(name+".Jet_SoftEl["+name+".nJet]/F").c_str());

    tree->Branch((name+".Jet_hist1").c_str(),       Jet_hist1       ,(name+".Jet_hist1["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_hist2").c_str(),       Jet_hist2       ,(name+".Jet_hist2["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_hist3").c_str(),       Jet_hist3       ,(name+".Jet_hist3["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_histJet").c_str(),     Jet_histJet     ,(name+".Jet_histJet["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_histSvx").c_str(),     Jet_histSvx     ,(name+".Jet_histSvx["+name+".nJet]/I").c_str());

    tree->Branch((name+".Jet_nFirstTrack").c_str(), Jet_nFirstTrack ,(name+".Jet_nFirstTrack["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_nLastTrack").c_str(),  Jet_nLastTrack  ,(name+".Jet_nLastTrack["+name+".nJet]/I").c_str()); 
    tree->Branch((name+".Jet_nFirstSV").c_str(),    Jet_nFirstSV    ,(name+".Jet_nFirstSV["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_nLastSV").c_str(),     Jet_nLastSV     ,(name+".Jet_nLastSV["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_SV_multi").c_str(),    Jet_SV_multi      ,(name+".Jet_SV_multi["+name+".nJet]/I").c_str());  
    tree->Branch((name+".Jet_nFirstTrkInc").c_str(), Jet_nFirstTrkInc ,(name+".Jet_nFirstTrkInc["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_nLastTrkInc").c_str(),  Jet_nLastTrkInc  ,(name+".Jet_nLastTrkInc["+name+".nJet]/I").c_str()); 

    tree->Branch((name+".Jet_VtxCat").c_str(),      Jet_VtxCat    ,(name+".Jet_VtxCat["+name+".nJet]/I").c_str() );
    tree->Branch((name+".Jet_looseID").c_str(),      Jet_looseID  ,(name+".Jet_looseID["+name+".nJet]/I").c_str());
    tree->Branch((name+".Jet_tightID").c_str(),      Jet_tightID  ,(name+".Jet_tightID["+name+".nJet]/I").c_str());
  } 

  void ReadTree(TTree *tree, std::string name="JetInfo") { 
    tree->SetBranchAddress((name+".Jet_pt").c_str(),          Jet_pt	       );
    tree->SetBranchAddress((name+".Jet_genpt").c_str(),       Jet_genpt       );
    tree->SetBranchAddress((name+".Jet_residual").c_str(),    Jet_residual    ); 
    tree->SetBranchAddress((name+".Jet_jes").c_str(),         Jet_jes         );
    tree->SetBranchAddress((name+".Jet_eta").c_str(),         Jet_eta         );
    tree->SetBranchAddress((name+".Jet_phi").c_str(),         Jet_phi         );
    tree->SetBranchAddress((name+".Jet_ntracks").c_str(),     Jet_ntracks     );
    tree->SetBranchAddress((name+".Jet_flavour").c_str(),     Jet_flavour     );
    //  tree->SetBranchAddress(name+".Jet_Ip1N",	   Jet_Ip1N	   );
    //  tree->SetBranchAddress(name+".Jet_Ip1P",	   Jet_Ip1P	   );
    tree->SetBranchAddress((name+".Jet_Ip2N").c_str(),        Jet_Ip2N        );
    tree->SetBranchAddress((name+".Jet_Ip2P").c_str(),        Jet_Ip2P        );
    tree->SetBranchAddress((name+".Jet_Ip3N").c_str(),        Jet_Ip3N        );
    tree->SetBranchAddress((name+".Jet_Ip3P").c_str(),        Jet_Ip3P        );
    //  tree->SetBranchAddress(name+".Jet_Ip4N",	   Jet_Ip4N	   );
    //  tree->SetBranchAddress(name+".Jet_Ip4P",	   Jet_Ip4P	   );
    //  tree->SetBranchAddress(name+".Jet_Mass4N",	   Jet_Mass4N	   );
    //  tree->SetBranchAddress(name+".Jet_Mass4P",	   Jet_Mass4P	   );
    tree->SetBranchAddress((name+".Jet_ProbaN").c_str(),      Jet_ProbaN      );
    tree->SetBranchAddress((name+".Jet_ProbaP").c_str(),      Jet_ProbaP      );
    //$$  tree->SetBranchAddress((name+".Jet_Proba").c_str(),       Jet_Proba       );
    tree->SetBranchAddress((name+".Jet_BprobN").c_str(),      Jet_BprobN      );
    tree->SetBranchAddress((name+".Jet_BprobP").c_str(),      Jet_BprobP      );
    //$$  tree->SetBranchAddress((name+".Jet_Bprob").c_str(),       Jet_Bprob       );
    tree->SetBranchAddress((name+".Jet_SvxN").c_str(),        Jet_SvxN        );
    tree->SetBranchAddress((name+".Jet_Svx").c_str(),         Jet_Svx         );
    //  tree->SetBranchAddress((name+".Jet_SvxNTracks").c_str(),  Jet_SvxNTracks  );
    //  tree->SetBranchAddress((name+".Jet_SvxTracks").c_str(),   Jet_SvxTracks   );
    tree->SetBranchAddress((name+".Jet_SvxNHP").c_str(),      Jet_SvxNHP      );
    tree->SetBranchAddress((name+".Jet_SvxHP").c_str(),       Jet_SvxHP       );
    tree->SetBranchAddress((name+".Jet_SvxMass").c_str(),     Jet_SvxMass     );
    tree->SetBranchAddress((name+".Jet_CombSvxN").c_str(),    Jet_CombSvxN    );
    tree->SetBranchAddress((name+".Jet_CombSvxP").c_str(),    Jet_CombSvxP    );
    tree->SetBranchAddress((name+".Jet_CombSvx").c_str(),     Jet_CombSvx     );

    tree->SetBranchAddress((name+".Jet_RetCombSvxN").c_str(), Jet_RetCombSvxN    );
    tree->SetBranchAddress((name+".Jet_RetCombSvxP").c_str(), Jet_RetCombSvxP    );
    tree->SetBranchAddress((name+".Jet_RetCombSvx").c_str(), Jet_RetCombSvx     );

    tree->SetBranchAddress((name+".Jet_CombCSVJP_N").c_str(), Jet_CombCSVJP_N    );
    tree->SetBranchAddress((name+".Jet_CombCSVJP_P").c_str(), Jet_CombCSVJP_P    );
    tree->SetBranchAddress((name+".Jet_CombCSVJP").c_str(), 	 Jet_CombCSVJP      );

    tree->SetBranchAddress((name+".Jet_CombCSVSL_N").c_str(), Jet_CombCSVSL_N    );
    tree->SetBranchAddress((name+".Jet_CombCSVSL_P").c_str(), Jet_CombCSVSL_P    );
    tree->SetBranchAddress((name+".Jet_CombCSVSL").c_str(), Jet_CombCSVSL      );

    tree->SetBranchAddress((name+".Jet_CombCSVJPSL_N").c_str(), Jet_CombCSVJPSL_N  );
    tree->SetBranchAddress((name+".Jet_CombCSVJPSL_P").c_str(), Jet_CombCSVJPSL_P  );
    tree->SetBranchAddress((name+".Jet_CombCSVJPSL").c_str(), Jet_CombCSVJPSL    );

    tree->SetBranchAddress((name+".Jet_SimpIVF_HP").c_str(),  Jet_SimpIVF_HP  );
    tree->SetBranchAddress((name+".Jet_SimpIVF_HE").c_str(),  Jet_SimpIVF_HE  );
    tree->SetBranchAddress((name+".Jet_DoubIVF_HE").c_str(),  Jet_DoubIVF_HE  );
    tree->SetBranchAddress((name+".Jet_CombIVF").c_str(),     Jet_CombIVF     );
    tree->SetBranchAddress((name+".Jet_CombIVF_P").c_str(), Jet_CombIVF_P   );

    tree->SetBranchAddress((name+".Jet_SoftMuN").c_str(), Jet_SoftMuN     );
    tree->SetBranchAddress((name+".Jet_SoftMuP").c_str(), Jet_SoftMuP     );
    tree->SetBranchAddress((name+".Jet_SoftMu").c_str(),  Jet_SoftMu      );

    tree->SetBranchAddress((name+".Jet_SoftElN").c_str(),     Jet_SoftElN     );
    tree->SetBranchAddress((name+".Jet_SoftElP").c_str(),     Jet_SoftElP     );
    tree->SetBranchAddress((name+".Jet_SoftEl").c_str(),      Jet_SoftEl      );

    tree->SetBranchAddress((name+".Jet_hist1").c_str(),       Jet_hist1       );
    tree->SetBranchAddress((name+".Jet_hist2").c_str(),       Jet_hist2       );
    tree->SetBranchAddress((name+".Jet_hist3").c_str(),       Jet_hist3       );
    tree->SetBranchAddress((name+".Jet_histJet").c_str(),     Jet_histJet     );
    tree->SetBranchAddress((name+".Jet_histSvx").c_str(),     Jet_histSvx     );

    tree->SetBranchAddress((name+".Jet_nFirstTrack").c_str(), Jet_nFirstTrack );
    tree->SetBranchAddress((name+".Jet_nLastTrack").c_str(),  Jet_nLastTrack  ); 
    tree->SetBranchAddress((name+".Jet_nFirstSV").c_str(),    Jet_nFirstSV    );
    tree->SetBranchAddress((name+".Jet_nLastSV").c_str(),     Jet_nLastSV     );
    tree->SetBranchAddress((name+".Jet_SV_multi").c_str(),    Jet_SV_multi      );  
    tree->SetBranchAddress((name+".Jet_nFirstTrkInc").c_str(),Jet_nFirstTrkInc );
    tree->SetBranchAddress((name+".Jet_nLastTrkInc").c_str(), Jet_nLastTrkInc  ); 

    tree->SetBranchAddress((name+".Jet_VtxCat").c_str(),      Jet_VtxCat  );
    tree->SetBranchAddress((name+".Jet_looseID").c_str(),     Jet_looseID);
    tree->SetBranchAddress((name+".Jet_tightID").c_str(),     Jet_tightID);
  } 

}; 

#endif 
