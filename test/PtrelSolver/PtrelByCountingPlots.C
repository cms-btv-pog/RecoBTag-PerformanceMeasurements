void PtrelByCountingPlots (const char * mfile, const char * mcfile)
{
    char name[256];

    gROOT->SetStyle("Plain");

    TFile * measurement = new TFile(mfile, "READ");
    TFile * mctruth = new TFile(mcfile, "READ");

    TCanvas * c1 = new TCanvas("c1", "PtrelSolver apply to TC");
    TImage *img = TImage::Create();

    c1->Divide(3,2);

    c1->cd(1);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_TCL", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_TCL");

    c1->cd(2);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_TCM", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_TCM");

    c1->cd(3);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_TCT", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_TCT");

    c1->cd(4);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_TCL", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_TCL");

    c1->cd(5);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_TCM", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_TCM");

    c1->cd(6);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_TCT", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_TCT");

    sprintf(name, "%s.TC.png", mfile);
    img->FromPad(c1);
    img->WriteImage(name);

    TCanvas * c2 = new TCanvas("c2", "PtrelSolver apply to JP");

    c2->Divide(3,2);

    c2->cd(1);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_JPL", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_JPL");

    c2->cd(2);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_JPM", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_JPM");

    c2->cd(3);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_JPT", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_JPT");

    c2->cd(4);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_JPL", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_JPL");

    c2->cd(5);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_JPM", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_JPM");

    c2->cd(6);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_JPT", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_JPT");

    sprintf(name, "%s.JP.png", mfile);
    img->FromPad(c2);
    img->WriteImage(name);

    TCanvas * c3 = new TCanvas("c3", "PtrelSolver apply to SV");

    c3->Divide(2,2);

    c3->cd(1);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_SVM", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_SVM");

    c3->cd(2);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_SVT", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_SVT");

    c3->cd(3);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_SVM", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_SVM");

    c3->cd(4);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_SVT", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_SVT");

    sprintf(name, "%s.SV.png", mfile);
    img->FromPad(c3);
    img->WriteImage(name);

    TCanvas * c4 = new TCanvas("c4", "PtrelSolver apply to CSV");

    c4->Divide(3,2);

    c4->cd(1);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_CSVL", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_CSVL");

    c4->cd(2);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_CSVM", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_CSVM");

    c4->cd(3);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_CSVT", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_CSVT");

    c4->cd(4);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_CSVL", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_CSVL");

    c4->cd(5);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_CSVM", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_CSVM");

    c4->cd(6);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_CSVT", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_CSVT");

    sprintf(name, "%s.CSV.png", mfile);
    img->FromPad(c4);
    img->WriteImage(name);

    TCanvas * c5 = new TCanvas("c5", "PtrelSolver apply to JBP");

    c5->Divide(3,2);

    c5->cd(1);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_JBPL", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_JBPL");

    c5->cd(2);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_JBPM", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_JBPM");

    c5->cd(3);
    plot(measurement, "/measurements/measurement_n_pT_ptag_pT_JBPT", mctruth, "/mctruth/mctruth_p_pT_b_ptag_pT_b_JBPT");

    c5->cd(4);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_JBPL", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_JBPL");

    c5->cd(5);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_JBPM", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_JBPM");

    c5->cd(6);
    plot(measurement, "/measurements/measurement_n_eta_ptag_eta_JBPT", mctruth, "/mctruth/mctruth_p_eta_b_ptag_eta_b_JBPT");

    sprintf(name, "%s.JBP.png", mfile);
    img->FromPad(c5);
    img->WriteImage(name);

}

void plot(TFile * measurement, const char * mname, TFile * mctruth, const char * mcname)
{
    TH1D * mHistogram = (TH1D*) measurement->Get(mname);
    TH1D * mcHistogram = (TH1D*) mctruth->Get(mcname);

    mHistogram->GetYaxis()->SetRangeUser(0., 1.2);
    mcHistogram->GetYaxis()->SetRangeUser(0., 1.2);

    mHistogram->SetLineWidth(2);
    mcHistogram->SetLineWidth(2);
    mcHistogram->SetLineStyle(2);
    mcHistogram->SetLineColor(kRed);

    mHistogram->Draw();
    mcHistogram->Draw("same");

    TLegend * legend = new TLegend(0.25, 0.77, 0.6, 0.92);
    legend->AddEntry(mHistogram, "Measure", "pl");
    legend->AddEntry(mcHistogram, "MCTruth", "l");
    legend->Draw();
}

