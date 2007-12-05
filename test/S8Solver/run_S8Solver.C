
void run_S8Solver() {

	//gROOT->SetStyle("CMS");
	gSystem->Load("libS8Solver");
	S8Solver s;
	s.SetData("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/histogramsV4/results_soup1_TCL_AwayTCL.root");
	s.SetName("TCL");
	s.Solve();
	s.Draw();
	
}
