
// PtrelSolver
// Author: Victor E. Bazterra, UIC (2008)

#include "PtrelBySystem4.h"

#include <math.h>
#include <iostream>

#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TFile.h"
#include "TPRegexp.h"

#include "../PtrelSolver/PtrelUtils.h"

ClassImp(PtrelBySystem4)

void PtrelBySystem4::solve(char const * inputfile, char const * outputfile)
{
    Info(__FUNCTION__, "Measuring efficiency");

    // Open a input file
    CreateSafely(TFile, input, TFile::Open(inputfile, "READ"))

    // Open a output file
    CreateSafely(TFile, output, TFile::Open(outputfile, "RECREATE"))

    // Creating sub directory for mesurements results
    if (!output->FindKey("measurements")) output->mkdir("measurements");

    char name[256];

    for (Int_t i = 1; i<Dependency::Dimension; ++i)
    {
        // Collection of flavor content n in sample
        ValueVector nValues;
        ValueVector nErrors;

        // Name for the n samples : n_(pT|eta)
        sprintf(name, "%s/n_%s", directory, Dependency::Name[i]);
        // Measure the flavor content in n sample
        CallSafely( measure(input, output, name, nValues, nErrors) )

        // Collection of flavor content p in sample
        ValueVector pValues;
        ValueVector pErrors;

        // Name for the p samples : p_(pT|eta)
        sprintf(name, "%s/p_%s", directory, Dependency::Name[i]);
        // Measure the flavor content in n sample
        CallSafely( measure(input, output, name, pValues, pErrors) )

        // Collection of counting in ntag samples
        StringVector ntagNames;
        ValueVector ntagValues;
        ValueVector ntagErrors;

        // Name for the ntag samples : ntag_(pT|eta)_[A-Z]*$
        sprintf(name, "ntag_%s_(?!SMT$)[A-Z]*$", Dependency::Name[i]);
        // Count the number of entries in ntag samples
        CallSafely( count(input, output, TPRegexp(name), ntagNames, ntagValues, ntagErrors) )

        // Collection of counting in ptag samples
        StringVector ptagNames;
        ValueVector ptagValues;
        ValueVector ptagErrors;

        // Name for the ptag samples : ptag_(pT|eta)_[A-Z]*$
        sprintf(name, "ptag_%s_(?!SMT$)[A-Z]*$", Dependency::Name[i]);
        // Count the number of entries in ptag samples
        CallSafely( count(input, output, TPRegexp(name), ptagNames, ptagValues, ptagErrors) )

        // HACK to force reading the objets from file
        // Close output
        input->Close();
        // Reopen
        GetSafely(input, TFile::Open(inputfile, "READ"))
        // ENDHACK

        // Get a 2D histogram to produce a place holder for the results later
        sprintf(name, "%s/n_%s", directory, Dependency::Name[i]);
        CreateSafely(TH2D, histogram2D, readTH2(input, name))

        for (std::size_t j = 0; j < ptagNames.size(); ++j)
        {
            // Collect the tag names
            TPRegexp t("[A-Z]{3,}");

            Int_t inx = ntagNames[j].Index(t);
            std::string tmp(ntagNames[j].Data());
            std::string tagA(tmp.substr(inx, tmp.size()));

            inx = ptagNames[j].Index(t);
            tmp = ptagNames[j].Data();
            std::string tagB(tmp.substr(inx, tmp.size()));

            // Read the beta coefficients
            sprintf(
                name, "coefficients/coefficient_mctruth_n_%s_b_ntag_%s_b_%s_mctruth_p_%s_b_ptag_%s_b_%s",
                Dependency::Name[i], Dependency::Name[i], tagA.c_str(), Dependency::Name[i], Dependency::Name[i], tagB.c_str()
            );

            std::cout << name << std::endl;

            CreateSafely( TH1D, betaHistogram, templates_->Get(name) )

            // Read the beta coefficients
            sprintf(
                name, "coefficients/coefficient_mctruth_n_%s_cl_ntag_%s_cl_%s_mctruth_p_%s_cl_ptag_%s_cl_%s",
                Dependency::Name[i], Dependency::Name[i], tagA.c_str(), Dependency::Name[i], Dependency::Name[i], tagB.c_str()
            );

            std::cout << name << std::endl;

            CreateSafely( TH1D, alfaHistogram, templates_->Get(name) )

            // Measure efficiencies histogram
            sprintf(name, "measurement_%s_%s_%s", ntagNames[j].Data(), ptagNames[j].Data(), Flavor::Name[Flavor::b]);
            Info(__FUNCTION__, "Computing efficiency %s", name);

            // This projection is just to get a place holder for the efficiencies
            TH1D * histogram1D = histogram2D->ProjectionX(name, -1, -1, "e");
            efficiencyHistogramSetup(histogram1D);

            CallSafely(
                compute(
                    histogram1D, nValues, nErrors,
                    pValues, pErrors,
                    ntagValues[j], ntagErrors[j],
                    ptagValues[j], ptagErrors[j],
                    alfaHistogram, betaHistogram
                )
            )

            // Saving the histogram
            output->cd();
            output->cd("measurements");
            histogram1D->Write();
        }
    }

    input->Close();
    output->Close();
}


bool PtrelBySystem4::compute(
    TH1 * histogram,
    ValueVector const & nValues,
    ValueVector const & nErrors,
    ValueVector const & pValues,
    ValueVector const & pErrors,
    TVectorD const & ntagValues,
    TVectorD const & ntagErrors,
    TVectorD const & ptagValues,
    TVectorD const & ptagErrors,
    TH1 * alfaHistogram,
    TH1 * betaHistogram
)
{
    // Loop over different bins
    for (std::size_t i = 0; i < nValues.size(); ++i)
    {
        Double_t nb = nValues[i](Flavor::b);
        Double_t nbError = nErrors[i](Flavor::b);
        Double_t n  = nValues[i].Sum();
        Double_t fnb = nb/n;
        Double_t fnbError = (n*nbError*nbError - nb*nb)/(n*n*n);

        Double_t pb = pValues[i](Flavor::b);
        Double_t pbError = pErrors[i](Flavor::b);
        Double_t p = pValues[i].Sum();
        Double_t fpb = pb/p;
        Double_t fpbError = (p*pbError*pbError - pb*pb)/(p*p*p);

        std::cout << "nb  : " << nb << " +- " << nbError << std::endl;
        std::cout << "pb  : " << pb << " +- " << pbError << std::endl;
        std::cout << "fnb : " << fnb << " +- " << sqrt(fnbError) << std::endl;
        std::cout << "fpb : " << fpb << " +- " << sqrt(fpbError) << std::endl;

        Double_t ntag = ntagValues(i);
        Double_t epsntag = ntag/n;
        Double_t epsntagError = epsntag*(1-epsntag)/n;

        Double_t ptag = ptagValues(i);
        Double_t epsptag = ptag/p;
        Double_t epsptagError = epsptag*(1-epsptag)/p;

        std::cout << "epsntag: " << epsntag << " +- " << sqrt(epsntagError) << std::endl;
        std::cout << "epsptag: " << epsptag << " +- " << sqrt(epsptagError) << std::endl;

        Double_t alpha = alfaHistogram->GetBinContent(i+1);
        Double_t alphaError = alfaHistogram->GetBinError(i+1);
        Double_t beta = betaHistogram->GetBinContent(i+1);
        Double_t betaError = alfaHistogram->GetBinError(i+1);

        std::cout << "alfa: " << alpha << " " << "beta: " << beta << std::endl;

        Double_t deltab = alpha*epsntag*(1.0-fpb) - epsptag*(1.0-fnb);
        Double_t delta  = alpha*fnb*(1.0-fpb) - beta*(1.0-fnb)*fpb;

        if ( delta == 0.0 )
        {
            histogram->SetBinContent(i+1, 0.);
            histogram->SetBinError(i+1, 0.);
            return true;
        }

        Double_t efficiency = deltab/delta;

        Double_t error = sqrt (
                             pow(epsptag/delta - deltab*(alpha*(1-fpb)+beta*fpb)/pow(delta, 2), 2) * fnbError +
                             pow(alpha*epsntag/delta - deltab*(alpha*fnb+beta*(1-fnb)/pow(delta, 2)), 2) * fpbError +
                             pow(alpha*(1-fpb)/delta, 2) * epsntagError +
                             pow((1-fnb)/delta, 2) * epsptagError +
                             pow((epsntag*(1-fpb)/delta - deltab*fnb*(1-fpb)/pow(delta,2)) * alphaError, 2) +
                             pow((deltab*(1-fnb)*fpb/pow(delta,2)) * betaError, 2)
                         );

        std::cout << "eff.: " << efficiency << " +- " << error << std::endl;

        histogram->SetBinContent(i+1, efficiency);
        histogram->SetBinError(i+1, error);
    }
    return true;
}


bool PtrelBySystem4::count(
    TFile * input,
    TFile * output,
    TPRegexp pattern,
    StringVector & hVector,
    ValueVector & vVector,
    ValueVector & eVector
) const
{
    // Return status
    bool status = false;

    // Clean the containers
    hVector.clear();
    vVector.clear();
    eVector.clear();

    // Move to the directory with 2d histograms
    input->cd(directory);

    // Loop over all keys in this directory
    TIter nextkey( gDirectory->GetListOfKeys() );

    TKey * key;

    char name[256];

    while (( key = (TKey*)nextkey() ))
    {
        // Measure those 2D histogram with name fit the pattern
        TObject * object = key->ReadObj();
        if ( object->IsA()->InheritsFrom( "TH2" ) )
            if ( TString(object->GetName() ).Contains(pattern) )
            {
                // Update status
                status = true;

                // Temporal container for efficiencies
                TH2 * histogram2D;

                // Process 2D histogram associated to the keyword
                GetSafelyZero(histogram2D, processTH2(object))

                // Get the number of bins in X direction
                Int_t nbins = histogram2D->GetNbinsX();

                // Temporal container with the results
                TVectorD values(nbins), errors(nbins);

                // Info
                Info(__FUNCTION__, "Counting the number of events in %s", histogram2D->GetName());

                // Loop over the bins
                for (Int_t i = 1; i <= nbins; ++i)
                {
                    // Set histogram name by the formula = data_x_bin (transient)
                    sprintf(name, "count_%s_%d", histogram2D->GetName(), i);
                    // Project histogram
                    TH1D * histogram1D = histogram2D->ProjectionY(name, i, i, "e");
                    // Appending the results
                    values(i-1) = histogram1D->Integral();
                    errors(i-1) = sqrt(values(i-1));
                }
                hVector.push_back(TString(histogram2D->GetName()));
                vVector.push_back(values);
                eVector.push_back(errors);
            }
    };

    if (!status) Error(__FUNCTION__, "Non matching histograms were found");

    return status;
}

