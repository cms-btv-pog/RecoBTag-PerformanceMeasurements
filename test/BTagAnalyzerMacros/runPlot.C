
{
gROOT->ProcessLine(".L DrawCommPlot.C++");

DrawMCgen("pt_hat","pt hat",1);
DrawMCgen("pt_hat_fin","pt hat",1);
DrawMCCompare("jet_ptgen","gen pT",1);
DrawMCgen("nPU_mc","nPU",1);
DrawMCgen("nPV_mc","nPV",1);

DrawMCCompare("jet_pt_all","reco pT",1);
DrawMCCompare("jet_eta","eta",0);

DrawMCCompare("CSV","CSV",1);
DrawMCCompare("TCHE","TCHE",1);
DrawMCCompare("JP","JP",1);
DrawMCCompare("CSVIVF","CSVv2 w/ IVF",1);

DrawMCCompare("track_multi"  ,      "number of tracks in the jets",0);
DrawMCCompare("trk_multi_sel"  ,    "number of selected tracks in the jets",0);


DrawMCgen("jet_pt_all_bfromg"   ,"pT of b jets from GSPL",1);

}
