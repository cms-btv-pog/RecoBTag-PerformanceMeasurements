
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

    char name[269], nsample[256];

    for (Int_t i = 1; i<Dependency::Dimension; ++i)
    {
        // Collection of flavor content n in sample
        ValueVector nValues;
        ValueVector nErrors;

        // Name for the n sample : n_(pT|eta)
        sprintf(name, "%s/n_%s", directory, Dependency::Name[i]);
        sprintf(nsample, "n_%s", Dependency::Name[i]);
        // Measure the flavor content in n sample
        CallSafely( measure(input, output, name, nValues, nErrors) )

        // Collection of flavor content ptag in samples
        ValueMatrix ptagValues;
        ValueMatrix ptagErrors;
        StringVector ptagHistograms;

        // Measure the flavor content in ptag samples
        sprintf(name, "ptag_%s_(?!SMT)[A-Z]*$", Dependency::Name[i]);
        CallSafely( measure(input, output, TPRegexp(name), ptagHistograms, ptagValues, ptagErrors) )

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
            CallSafely( compute(histogram1D, mistag1D, pValues, nValues, nErrors, ptagValues[j], ptagErrors[j]) )

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
    ValueVector const & nErrors,
    ValueVector const & ptagValues,
    ValueVector const & ptagErrors
)
{
    // Loop over different bins
    for (std::size_t i = 0; i < nValues.size(); ++i)
    {
        Double_t ptagb = ptagValues[i](Flavor::b);
        Double_t ptagbError = ptagErrors[i](Flavor::b);

        Double_t nb = nValues[i](Flavor::b);
        Double_t nbError = nErrors[i](Flavor::b);

        Double_t ncl = nValues[i].Sum() - nb;
        Double_t nclError = sqrt( nValues[i].Sum() + nbError * nbError );

        Double_t m = mistag->GetBinContent(i+1);
        Double_t p = pValues(i);

        Double_t efficiency, error;

        if (ncl != 0. && ptagb != 0.)
        {
            efficiency = ptagb/(p - m * ncl);
            error = sqrt(pow(ptagbError,2)/pow(p - m * ncl, 2) + pow(nclError,2)*pow(ptagb*m/pow(p - m * ncl, 2), 2));
        }
        else
        {
            efficiency = 0.;
            error = 0.;
        }

        histogram->SetBinContent(i+1, efficiency);
        histogram->SetBinError(i+1, error);
    }
    return true;
}
