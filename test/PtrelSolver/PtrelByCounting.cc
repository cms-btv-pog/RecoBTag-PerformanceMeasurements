
// PtrelSolver
// Author: Victor E. Bazterra, UIC (2008)

#include "PtrelByCounting.h"

#include <math.h>
#include <iostream>

#include "TFile.h"
#include "TPRegexp.h"

#include "PtrelUtils.h"

ClassImp(PtrelByCounting)

void PtrelByCounting::solve(char const * inputfile, char const * mistagfile, char const * outputfile)
{
    Info(__FUNCTION__, "Measuring efficiency");

    // Open a input file
    CreateSafely(TFile, input, TFile::Open(inputfile, "READ"))

    // Open the mistag rate
    CreateSafely(TFile, mistagf, TFile::Open(mistagfile, "READ"))

    // Open a output file
    CreateSafely(TFile, output, TFile::Open(outputfile, "RECREATE"))

    // Creating sub directory for templates
    if (!output->FindKey("measurements")) output->mkdir("measurements");

    char name[256], nsample[256];

    for (Int_t i = 1; i<Dependency::Dimension; ++i)
    {
        // Collection of flavor content n in sample
        ValueVector nValues;
        CovarianceVector nCovariances;

        // Name for the n sample : n_(pT|eta)
        sprintf(name, "%s/n_%s", directory, Dependency::Name[i]);
        sprintf(nsample, "n_%s", Dependency::Name[i]);
        // Measure the flavor content in n sample
        CallSafely( measure(input, output, name, nValues, nCovariances) )

        // Collection of flavor content ptag in samples
        ValueMatrix ptagValues;
        CovarianceMatrix ptagCovariances;
        StringVector ptagHistograms;

        // Measure the flavor content in ptag samples
        sprintf(name, "ptag_%s_[A-Z]*$", Dependency::Name[i]);
        CallSafely( measure(input, output, TPRegexp(name), ptagHistograms, ptagValues, ptagCovariances) )

        // HACK to force reading the objets from file
        // Close output
        input->Close();
        // Reopen
        GetSafely(input, TFile::Open(inputfile, "READ"))
        // ENDHACK

        // Get a 2D histogram to produce a place holder for the results later
        sprintf(name, "%s/n_%s", directory, Dependency::Name[i]);
        CreateSafely(TH2D, histogram2D, readTH2(input, name))

        // Total number of event in sample p (this information will be measure)
        sprintf(name, "%s/p_%s", directory, Dependency::Name[i]);
        CreateSafely(TH2D, pHistogram2D, readTH2(input, name))
        TVectorD pValues(pHistogram2D->GetNbinsX());

        for (Int_t j = 0; j < pHistogram2D->GetNbinsX(); ++j)
        {
            TH1 * pHistogram1D = pHistogram2D->ProjectionY(name, j+1, j+1, "e");
            pValues(j) = pHistogram1D->Integral();
        }

        // HACK: Get the mistag rate from file
        sprintf(name, "/mctruth/mctruth_n_%s_cl_p_%s_cl", Dependency::Name[i], Dependency::Name[i]);
        Info(__FUNCTION__, "Reading mistag rate info %s", name);
        CreateSafely(TH1D, mistag1D, mistagf->Get(name))

        for (std::size_t j = 0; j < ptagHistograms.size(); ++j)
        {
            // Measure efficiencies histogram
            sprintf(name, "measurement_%s_%s", nsample, ptagHistograms[j].Data());

            Info(__FUNCTION__, "Computing efficiency %s", name);

            // This projection is just to get a place holder for the efficiencies
            TH1D * histogram1D = histogram2D->ProjectionX(name, -1, -1, "e");
            efficiencyHistogramSetup(histogram1D);

            // Calculate b-efficiencies
            CallSafely( compute(histogram1D, mistag1D, pValues, nValues, nCovariances, ptagValues[j], ptagCovariances[j]) )

            // Saving the histogram
            output->cd();
            output->cd("measurements");
            histogram1D->Write();
        }
    }

    input->Close();
    output->Close();
}


bool PtrelByCounting::compute(
    TH1 * histogram,
    TH1 * mistag,
    TVectorD const & pValues,
    ValueVector const & nValues,
    CovarianceVector const & nCovariance,
    ValueVector const & ptagValues,
    CovarianceVector const & ptagCovariance
)
{
    // Loop over different bins
    for (std::size_t i = 0; i < nValues.size(); ++i)
    {
        Double_t n_cl = nValues[i](Flavor::cl);
        Double_t ptag_b = ptagValues[i](Flavor::b);

        Double_t n_cl_error = nCovariance[i](Flavor::cl, Flavor::cl);
        Double_t ptag_b_error = ptagCovariance[i](Flavor::b, Flavor::b);

        Double_t m = mistag->GetBinContent(i+1);
        Double_t p = pValues(i);

        histogram->SetBinContent(i+1, ptag_b/(p - m * n_cl));
        histogram->SetBinError(i+1,
                               sqrt(
                                   ptag_b_error/pow(p - m * n_cl, 2) + n_cl_error*pow(ptag_b*m/pow(p - m * n_cl, 2), 2)
                               )
                              );
    }
    return true;
}

