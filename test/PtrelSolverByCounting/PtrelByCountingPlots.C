void PtrelByCountingPlots (const char * mfile, const char * mcfile)
{
    gROOT->SetStyle("Plain");

    TFile * measurement = new TFile(mfile, "READ");
    TFile * mctruth = new TFile(mcfile, "READ");

    TCanvas * c1 = new TCanvas("c1", "PtrelByCounting apply to TC");

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


    TCanvas * c2 = new TCanvas("c2", "PtrelByFitting apply to JP");

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

