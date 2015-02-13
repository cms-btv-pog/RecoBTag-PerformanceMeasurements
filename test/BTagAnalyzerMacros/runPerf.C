{

gROOT->ProcessLine(".L ComparePerformance.C++");

DrawAllPlot(1,2);
DrawAllPlot(1,1);

DrawPlot("TCHE_1"            ,1,"QCD GEN 50-200 : TCHE");
DrawPlot("TCHP_1"            ,1,"QCD GEN 50-200 : TCHP");
DrawPlot("JP_1"              ,1,"QCD GEN 50-200 : JP");
DrawPlot("JBP_1"             ,1,"QCD GEN 50-200 : JBP");
DrawPlot("SSV_1"             ,1,"QCD GEN 50-200 : SSVHE");
DrawPlot("SSVHP_1"           ,1,"QCD GEN 50-200 : SSVHP");
DrawPlot("CSV_1"             ,1,"QCD GEN 50-200 : CSV");
DrawPlot("CSVIVF_1"          ,1,"QCD GEN 50-200 : CSVv2 w/ IVF");

DrawPlot("TCHE_2"            ,1,"QCD GEN 500-1000 : TCHE");
DrawPlot("TCHP_2"            ,1,"QCD GEN 500-1000 : TCHP");
DrawPlot("JP_2"              ,1,"QCD GEN 500-1000 : JP");
DrawPlot("JBP_2"             ,1,"QCD GEN 500-1000 : JBP");
DrawPlot("SSV_2"             ,1,"QCD GEN 500-1000 : SSVHE");
DrawPlot("SSVHP_2"           ,1,"QCD GEN 500-1000 : SSVHP");
DrawPlot("CSV_2"             ,1,"QCD GEN 500-1000 : CSV");
DrawPlot("CSVIVF_2"          ,1,"QCD GEN 500-1000 : CSVv2 w/ IVF");

DrawEff("CSV_1"             ,1,"QCD GEN 50-200 : CSV");
DrawEff("CSVIVF_1"             ,1,"QCD GEN 50-200 : CSVIVF");
DrawEff("TCHE_1"             ,1,"QCD GEN 50-200 : TCHE");
DrawEff("TCHP_1"             ,1,"QCD GEN 50-200 : TCHP");
DrawEff("JP_1"             ,1,"QCD GEN 50-200 : JP");
DrawEff("SSV_1"             ,1,"QCD GEN 50-200 : SSVHE");

DrawEffvsPt("tchpt", 1, "TCHP T vs pT", 1);
DrawEffvsPt("jpl", 1, "JP L vs pT", 1);
DrawEffvsPt("csvivfl", 1, "CSVv2 w/ IVF L vs pT", 1);
DrawEffvsPt("csvivfm", 1, "CSVv2 w/ IVF M vs pT", 1);
DrawEffvsPt("csvl", 0, "CSV L vs pT", 2);
DrawEffvsPt("csvl", 1, "CSV L vs pT", 2);
DrawEffvsPt("csvm", 1, "CSV M vs pT", 2);
DrawEffvsPt("csvivfl", 1, "CSVv2 w/ IVF L vs pT", 2);
DrawEffvsPt("tchpt", 1, "TCHP T vs pT", 2);
DrawEffvsPt("csvivfl", 1, "CSVv2 w/ IVF L vs pT", 2);
DrawEffvsPt("csvm", 1, "CSV M vs pT", 2);
DrawEffvsPt("csvl", 1, "CSV L vs pT", 2);

gROOT->ProcessLine(".L ComparePerformance.C++");
DrawEff("TCHP_1"             ,1,"QCD GEN 50-200 : TCHP");


}
