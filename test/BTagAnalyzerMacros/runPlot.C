
{
gROOT->ProcessLine(".L DrawCommPlot.C++");

DrawMCCgen("pt_hat","pt hat",1)
DrawMCgen("pt_hat","pt hat",1)
DrawMCgen("pt_hat_fin1","pt hat",1)
DrawMCgen("pt_hat_fin2","pt hat",1)
DrawMCCompare("jet1_ptgen","gen pT",1)
DrawMCCompare("jet2_ptgen","gen pT",1)
DrawMCgen("nPU_mc","nPU",1)
DrawMCgen("nPV_mc","nPV",1)

DrawMCCompare("jet2_pt_all","reco pT",1)
DrawMCCompare("jet1_pt_all","reco pT",1)
DrawMCCompare("jet2_eta","eta",1)
DrawMCCompare("jet2_eta","eta",0)
DrawMCCompare("jet1_eta","eta",0)
DrawMCCompare("jet1_diffpt","diff pt",0)
DrawMCCompare("jet1_diffrel","rel diff pt",0)
DrawMCCompare("jet1_diffrel","rel diff pt",1)
DrawMCCompare("jet2_diffrel","rel diff pt",1)
DrawMCCompare("jet2_diffpt","diff pt",1)
DrawMCCompare("CSV_1","CSV",1)
DrawMCCompare("TCHE_1","TCHE",1)
DrawMCCompare("JP_1","JP",1)
DrawMCCompare("CSVIVF_1","CSVv2 w/ IVF",1)


/*
gROOT->ProcessLine(".L DrawCommPlot.C++");
DrawMCgen("jet_pt_all_bfromg"   ,"pT of b jets from GSPL",1);
DrawMCgen("jet_pt_all_c"   ,"pT of c jets",1);
DrawMCgen("jet_pt_all_cfromg"   ,"pT of c jets from GSPL",1);
DrawMCgen("jet_pt_all_l"   ,"pT of l jets",1);
DrawMCgen("jet_pt_all_b"   ,"pT of b jets",1);
DrawMCgen("jet_pt_all_bfromg"   ,"pT of b jets from GSPL",1);
DrawMCgen("jet_pt_all_l"   ,"pT of l jets",1);
DrawMCgen("jet_pt_all_c"   ,"pT of c jets",1);
DrawMCgen("jet_pt_all_cfromg"   ,"pT of c jets from GSPL",1);
DrawMCCompare("jet_pt_all"   ,"pT of all jets",1);
DrawMCCompare("jet_pt_all"   ,"pT of all jets",0);
DrawMCgen("jet_pt_all_cfromg"   ,"pT of c jets from GSPL",0);
DrawMCgen("jet_pt_all_c"   ,"pT of c jets",0);
DrawMCgen("jet_pt_all_l"   ,"pT of l jets",0);
DrawMCgen("jet_pt_all_bfromg"   ,"pT of b jets from GSPL",0);
DrawMCgen("jet_pt_all_b"   ,"pT of b jets",0);
        DrawMCCompare("track_multi"  ,      "number of tracks in the jets",0);
        DrawMCCompare("trk_multi_sel"  ,    "number of selected tracks in the jets",0);
        DrawMCCompare("muon_multi"   ,      "number of muons", 1);
        DrawMCCompare("muon_multi_sel"   ,  "number of selected muons",1);
        DrawMCCompare("mu_ptrel"     ,      "p_{T} rel. of the muon",0);
        DrawMCCompare("mu_chi2"      ,      "norm. #chi^{2} of the muon", 1);
        DrawMCCompare("muon_Pt",             "Muon p_{T}",1);
        DrawMCCompare("muon_eta",            "Muon #eta",0);
        DrawMCCompare("muon_phi",            "Muon #phi",0);
        DrawMCCompare("muon_Ip3d",           "Muon 3D IP",1);
        DrawMCCompare("muon_Ip2d",           "Muon 2D IP",1);
        DrawMCCompare("muon_Sip3d",          "Muon 3D IP significance",1);
        DrawMCCompare("muon_Sip2d",          "Muon 2D IP significance",1);
        DrawMCCompare("muon_DeltaR",         "Muon #Delta R",0);
*/



}
