#include "PtrelByFitting.h"

#include <math.h>

#include "TFile.h"
#include "TPRegexp.h"

#include "PtrelUtils.h"

ClassImp(PtrelByFitting)

void PtrelByFitting::solve(char const * inputfile, char const * outputfile)
{
    Info(__FUNCTION__, "Measuring efficiency");

    // Opening a input file
    CreateSafely(TFile, input, TFile::Open(inputfile, "READ"))

    // Opening a output file
    CreateSafely(TFile, output, TFile::Open(outputfile, "RECREATE"))

    // Creating sub directory for templates
    if (!output->FindKey("measurements")) output->mkdir("measurements");

    char name[256];

    // Methods based on n and p sample
    char const * methods[] = {"n","p"};
    // Looping over the different methods
    for (Int_t m = 0; m < 2; ++m)
    {
        // Measuring over all the dependencies
        for (Int_t i = 1; i < Dependency::Dimension; ++i)
        {
            // Collection of flavor content mesurements for making the denominators
            StringVector dHistograms;
            ValueMatrix dValues;
            ValueMatrix dErrors;

            sprintf(name, "%s_%s", methods[m], Dependency::Name[i]);
            CallSafely( measure(input, output, TPRegexp(name), dHistograms, dValues, dErrors) )

            // Collection of flavor content mesurements for making the numerators
            StringVector nHistograms;
            ValueMatrix nValues;
            ValueMatrix nErrors;

            sprintf(name, "%stag_%s", methods[m], Dependency::Name[i]);
            CallSafely( measure(input, output, TPRegexp(name), nHistograms, nValues, nErrors) )

            // HACK to force reading the objets from file
            // Close input
            input->Close();
            // Reopen
            GetSafely(input, TFile::Open(inputfile, "READ"))
            // ENDHACK

            for (std::size_t j = 0; j < dHistograms.size(); ++j)
            {
                // Process the histogram rebinning etc.
                sprintf(name, "%s/%s", directory, dHistograms[j].Data());
                CreateSafely(TH2D, histogram2D, readTH2(input, name))

                for (std::size_t k = 0; k < nHistograms.size(); ++k)
                {
                    // Measure efficiencies histogram
                    sprintf(name, "measurement_%s_%s_%s", dHistograms[j].Data(), nHistograms[k].Data(), Flavor::Name[Flavor::b]);
                    Info(__FUNCTION__, "Computing efficiency %s", name);
                    // This projection is just to get a place holder for the efficiencies
                    TH1D * histogram1D = histogram2D->ProjectionX(name, -1, -1, "e");
                    efficiencyHistogramSetup(histogram1D);
                    // Calculate b-efficiencies
                    CallSafely( compute(histogram1D, Flavor::b, dValues[j], dErrors[j], nValues[k], nErrors[k]) )
                    // Saving the histogram
                    output->cd();
                    output->cd("measurements");
                    histogram1D->Write();
                }
            }
        }
    }

    input->Close();
    output->Close();
}

bool PtrelByFitting::compute(
    TH1 * histogram,
    Flavor::Type flavor,
    ValueVector const & dValues,
    ValueVector const & dErrors,
    ValueVector const & nValues,
    ValueVector const & nErrors
)
{
    // Loop over different bins
    for (std::size_t i = 0; i < dValues.size(); ++i)
    {
        Double_t d = dValues[i](flavor);
        Double_t n = nValues[i](flavor);

        TVectorD errors(2);
        errors(0) = dErrors[i](flavor);
        errors(1) = nErrors[i](flavor);

        TVectorD derivatives(2);
        derivatives(0) = - n/(d*d);
        derivatives(1) = 1./n;

        kernel_.removeCorrelations();

        if (d!=0)
        {
            histogram->SetBinContent(i+1, n/d);
            histogram->SetBinError(i+1, kernel_(derivatives, errors));
        }
        else
        {
            histogram->SetBinContent(i+1,0);
            histogram->SetBinError(i+1,0);
        }
    }
    return true;
}

