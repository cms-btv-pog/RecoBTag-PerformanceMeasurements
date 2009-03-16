void DumpPtrelByFittingMeasurements (const char * mfile, const char * tfile)
{
    char name[256];

    gROOT->SetStyle("Plain");

    TFile * measurement = new TFile(mfile, "READ");

    std::ostringstream tables;

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_TCL_b",
        "/measurements/measurement_n_eta_ntag_eta_TCL_b",
        "TrackCountingHighEff_Loose",
        2.0
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_TCM_b",
        "/measurements/measurement_n_eta_ntag_eta_TCM_b",
        "TrackCountingHighEff_Medium",
        4.2
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_TCT_b",
        "/measurements/measurement_n_eta_ntag_eta_TCT_b",
        "TrackCountingHighPur_Tight",
        4.1
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_JPL_b",
        "/measurements/measurement_n_eta_ntag_eta_JPL_b",
        "JetProbability_Loose",
        0.24
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_JPM_b",
        "/measurements/measurement_n_eta_ntag_eta_JPM_b",
        "JetProbability_Medium",
        0.49
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_JPT_b",
        "/measurements/measurement_n_eta_ntag_eta_JPT_b",
        "JetProbability_Tight",
        0.74
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_SVM_b",
        "/measurements/measurement_n_eta_ntag_eta_SVM_b",
        "SimpleSecondaryVertex_Medium",
        2.1
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_SVT_b",
        "/measurements/measurement_n_eta_ntag_eta_SVT_b",
        "SimpleSecondaryVertex_Tight",
        3.6
    );


    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_CSVL_b",
        "/measurements/measurement_n_eta_ntag_eta_CSVL_b",
        "CombinedSecondaryVertex_Loose",
        0.39
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_CSVM_b",
        "/measurements/measurement_n_eta_ntag_eta_CSVM_b",
        "CombinedSecondaryVertex_Medium",
        0.84
    );

    PtrelMeasurementTable(
        tables,
        measurement,
        "/measurements/measurement_n_pT_ntag_pT_CSVT_b",
        "/measurements/measurement_n_eta_ntag_eta_CSVT_b",
        "CombinedSecondaryVertex_Tight",
        0.95
    );

    ofstream output(tfile);
    output << tables.str();
}

void PtrelMeasurementTable(std::ostringstream & results, TFile * measurement, const char * ptdep, const char * etadep, const char * tagger, Double_t cut)
{
    TH1D * ptHistogram = (TH1D*) measurement->Get(ptdep);
    TH1D * etaHistogram = (TH1D*) measurement->Get(etadep);

    // Output stream with the table
    std::ostringstream results;

    // Header of for the table
    results << tagger << endl << cut << endl;
    results << "BtagPerformancePayloadFromTableEtaJetEt" << endl << '6' << endl;
    results << "etamin_etamax_etmin_etmax_bef_berr" << endl;

    // Get apointer to eta and pt axis
    TAxis * etaAxis = etaHistogram->GetXaxis();
    TAxis * ptAxis = ptHistogram->GetXaxis();

    Double_t eff_pt, eff_eta, eff;
    Double_t eff_pt_error, eff_eta_error, eff_error;

    for (Int_t etaIndex = 1; etaIndex <= etaAxis->GetNbins(); ++etaIndex)
        for (Int_t ptIndex = 1; ptIndex <= ptAxis->GetNbins(); ++ptIndex)
        {
            // Axis values
            results << etaAxis->GetBinLowEdge(etaIndex) << " ";
            results << etaAxis->GetBinUpEdge(etaIndex) << " ";
            results << ptAxis->GetBinLowEdge(ptIndex) << " ";
            results << ptAxis->GetBinUpEdge(ptIndex) << " ";

            // Getting the efficiencies
            eff_eta = etaHistogram->GetBinContent(etaIndex);
            eff_pt = ptHistogram->GetBinContent(ptIndex);

            eff = eff_eta*eff_pt;

            eff_error = 0;

            // Error contributed by eta dependency
            eff_eta_error = ptHistogram->GetBinError(etaIndex);
            eff_error += pow(eff_pt, 2) * pow(eff_eta_error, 2);

            // Error contributed by pt dependency
            eff_pt_error = ptHistogram->GetBinError(ptIndex);
            eff_error += pow(eff_eta, 2) * pow(eff_pt_error, 2);

            // Total error
            eff_error = sqrt(eff_error);

            // Streaming the results
            results << eff << " " << eff_error << endl;
        }
    results << endl;
}

