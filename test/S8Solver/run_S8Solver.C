
void run_S8Solver() {

	//gROOT->SetStyle("CMS");
	gSystem->Load("libS8Solver");
	S8Solver s;
	s.SetData("mc_ssvhp_t_31aug_pf.root");
	s.SetCorrData("mc_ssvhp_t_31aug_pf.root");
	s.SetName("TCHEL");
	s.SetAlphaConstant(false);
	s.SetBetaConstant(true);
	s.SetGammaConstant(true);
	s.SetDeltaConstant(false);
	s.SetKappabConstant(true);
	s.SetKappaclConstant(false);
	s.UseMCTrue(false);

	s.Solve();
	s.Draw();
	s.PrintData();
	s.PrintData("input");
	
//	s.Draw();
	
}

