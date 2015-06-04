
{
gROOT->ProcessLine(".L DrawCommPlot.C++");

DrawMCgen("pt_hat","pt hat",1);
DrawMCgen("pt_hat_fin","pt hat",1);
DrawMCgen("nPU_mc","nPU",1);
DrawMCgen("nPV_mc","nPV",1);
DrawMCgen("jet_pt_all_bfromg"   ,"pT of b jets from GSPL",1);

DrawMC("jet_ptgen","gen pT",1,89731000000);
DrawMC("jet_pt_all","reco pT",1,89731000000);
DrawMC("jet_eta","eta",0,89731000000);

DrawMC("CSV","CSV",1,89731000000);
DrawMC("TCHE","TCHE",1,89731000000);
DrawMC("JP","JP",1,89731000000);
DrawMC("CSVIVF","CSVv2 w/ IVF",1,89731000000);

DrawMC("track_multi"  ,      "number of tracks in the jets",0,89731000000);
//DrawMCCompare("trk_multi_sel"  ,    "number of selected tracks in the jets",0);

}
