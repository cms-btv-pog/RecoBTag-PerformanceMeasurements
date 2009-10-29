
#include "PtrelTemplateMaker.h"

#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjString.h"

#include "PtrelUtils.h"

ClassImp(PtrelTemplateMaker)


char const * PtrelTemplateMaker::directory = "/MCTruth";


void PtrelTemplateMaker::make(
    char const * inputfile,
    char const * outputfile
) const
{
    // Information
    Info(__FUNCTION__, "Reading %s input file", inputfile);

    // Opening a input file
    CreateSafely(TFile, input, TFile::Open(inputfile, "READ"))

    // Opening a output file
    CreateSafely(TFile, output, TFile::Open(outputfile, "RECREATE"))

    // Make templates
    CallSafely( makeTemplates(input, output) )

    char dName[256], nName[256];

    for (long i = 1; i < Dependency::Dimension; ++i)
    {
        // Make n-effciencies
        sprintf(dName, "n_%s_b", Dependency::Name[i]);
        sprintf(nName, "ntag_%s_b$", Dependency::Name[i]);
        CallSafely( makeEfficiencies(input, output, TPRegexp(dName), TPRegexp(nName)) )

        // Make p-efficiencies
        sprintf(dName, "p_%s_b", Dependency::Name[i]);
        sprintf(nName, "ptag_%s_b$", Dependency::Name[i]);
        CallSafely( makeEfficiencies(input, output, TPRegexp(dName), TPRegexp(nName)) )

        // Make n-cl-effciencies
        sprintf(dName, "n_%s_cl", Dependency::Name[i]);
        sprintf(nName, "ntag_%s_cl", Dependency::Name[i]);
        CallSafely( makeEfficiencies(input, output, TPRegexp(dName), TPRegexp(nName)) )

        // Make p-cl-efficiencies
        sprintf(dName, "p_%s_cl", Dependency::Name[i]);
        sprintf(nName, "ptag_%s_cl", Dependency::Name[i]);
        CallSafely( makeEfficiencies(input, output, TPRegexp(dName), TPRegexp(nName)) )

        // Make mistag for counting
        sprintf(dName, "n_%s_cl", Dependency::Name[i]);
        sprintf(nName, "p_%s_cl", Dependency::Name[i]);
        CallSafely( makeEfficiencies(input, output, TPRegexp(dName), TPRegexp(nName)) )

        // Evaluate the coefficients between b-efficiencies from p and n samples
        sprintf(dName, "mctruth_n_%s_b_ntag_%s_b", Dependency::Name[i], Dependency::Name[i]);
        sprintf(nName, "mctruth_p_%s_b_ptag_%s_b", Dependency::Name[i], Dependency::Name[i]);
        CallSafely( makeCoefficients(output, TPRegexp(dName), TPRegexp(nName)) )

        // Evaluate the coefficients between cl-efficiencies from p and n samples
        sprintf(dName, "mctruth_n_%s_cl_ntag_%s", Dependency::Name[i], Dependency::Name[i]);
        sprintf(nName, "mctruth_p_%s_cl_ptag_%s", Dependency::Name[i], Dependency::Name[i]);
        CallSafely( makeCoefficients(output, TPRegexp(dName), TPRegexp(nName)) )
    }

    // Closing the files
    input->Close();
    output->Close();
}


TH2* PtrelTemplateMaker::processTH2(TObject * object) const
{
    // Cast the object pointer into 2D histogram
    TH2 * histogram2D = (TH2*) object;

    // Check for ptrel rebinning
    if ( rebin_[Dependency::ptrel] != 1 )
    {
        Info(__FUNCTION__, "Rebinning ptrel a factor %d in %s", rebin_[Dependency::ptrel], histogram2D->GetName());
        histogram2D->RebinY( rebin_[Dependency::ptrel] );
    }

    // Check the histogram dependency and rebin
    for (Int_t j = 1; j < Dependency::Dimension; ++j)
        if ( TString(histogram2D->GetName()).Contains(Dependency::Name[j]) && rebin_[j] != 1)
        {
            Info(__FUNCTION__, "Rebinning %s a factor %d in %s", Dependency::Name[j], rebin_[j], histogram2D->GetName());
            histogram2D->RebinX( rebin_[j] );
        }

    return histogram2D;
}


bool PtrelTemplateMaker::makeCoefficients (
    TFile * output,
    TPRegexp patternD,
    TPRegexp patternN
) const
{
    char name[256];

    // Return status
    bool status = false;

    // Information
    Info(__FUNCTION__, "Starting making coefficients");

    // Creating sub directory for coefficients between mctruth efficiencies
    if (!output->FindKey("coefficients")) output->mkdir("coefficients");

    // Move to the directory with 2d histograms
    output->cd("/mctruth");

    // Loop over all denominator keys in this directory
    TIter nextkeyD( gDirectory->GetListOfKeys() );

    // Loop over all numerator keys in this directory
    TIter nextkeyN( gDirectory->GetListOfKeys() );

    // Iterator for the denominators
    TKey * keyD;
    TKey * keyN;

    while (( keyD = (TKey*)nextkeyD() ))
    {
        // Select only 2D histograms
        TObject * objectD = keyD->ReadObj();
        if ( objectD->IsA()->InheritsFrom( "TH1" ) )
            // Select those histogram that match the pattern
            if ( TString(objectD->GetName()).Contains(patternD) )
            {
                // Information
                Info(__FUNCTION__, "Selecting as denominator %s", objectD->GetName());

                // Cast the object pointer into 2D histogram
                TH1D * denominator = (TH1D*) objectD;

                while (( keyN = (TKey*)nextkeyN() ))
                {
                    // Select only 2D histograms
                    TObject * objectN = keyN->ReadObj();
                    if ( objectN->IsA()->InheritsFrom( "TH1" ) )
                        // Select those histogram that match the pattern
                        if ( TString(objectN->GetName()).Contains(patternN) )
                        {
                            // Update status
                            status = true;

                            // Information
                            Info(__FUNCTION__, "Selecting as numerator %s", objectN->GetName());

                            // Cast the object pointer into 2D histogram
                            TH1D * numerator = (TH1D*) objectN;

                            // MCTruth efficiencies histogram
                            sprintf(name, "coefficient_%s_%s", denominator->GetName(), numerator->GetName());
                            Info(__FUNCTION__, "Calculating coefficient %s", name);
                            TH1D * coefficient = (TH1D*) numerator->Clone();
                            coefficient->SetName(name);
                            coefficient->Divide(numerator, denominator, 1., 1., "e");
                            efficiencyHistogramSetup(coefficient);

                            // Save the beta coefficient histogram
                            output->cd("coefficients");
                            coefficient->Write();
                            break;
                        }
                };
            }
    };

    if (!status) Error(__FUNCTION__, "Non matching histograms were found");

    return status;
}


bool PtrelTemplateMaker::makeEfficiencies (
    TFile * input,
    TFile * output,
    TPRegexp patternD,
    TPRegexp patternN
) const
{
    char name[256];

    // Return value
    bool status = false;

    // Information
    Info(__FUNCTION__, "Starting making efficiencies");

    // Creating sub directory for templates
    if (!output->FindKey("mctruth")) output->mkdir("mctruth");

    // Move to the directory with 2d histograms
    input->cd(PtrelTemplateMaker::directory);

    // Loop over all denominator keys in this directory
    TIter nextkeyD( gDirectory->GetListOfKeys() );

    // Iterator for the denominators
    TKey * keyD;

    while (( keyD = (TKey*)nextkeyD() ))
    {
        // Select only 2D histograms
        TObject * objectD = keyD->ReadObj();
        if ( objectD->IsA()->InheritsFrom( "TH2" ) )
            // Select those histogram that match the pattern
            if ( TString(objectD->GetName()).Contains(patternD) )
            {
                // Updating status
                status = true;

                // Information
                Info(__FUNCTION__, "Selecting as denominator %s", objectD->GetName());

                // Cast the object pointer into 2D histogram
                TH2D * denominator2D = (TH2D*) processTH2(objectD);

                // Denominator histogram
                sprintf(name, "denominator_%s", denominator2D->GetName());
                TH1D * denominator1D = denominator2D->ProjectionX(name, -1, -1, "e");

                // Loop over all numerator keys in this directory
                TIter nextkeyN( gDirectory->GetListOfKeys() );

                // Iterator for the numerator
                TKey * keyN;

                while (( keyN = (TKey*)nextkeyN() ))
                {
                    // Select only 2D histograms
                    TObject * objectN = keyN->ReadObj();
                    if ( objectN->IsA()->InheritsFrom( "TH2" ) )
                        // Select those histogram that match the pattern
                        if ( TString(objectN->GetName()).Contains(patternN) )
                        {
                            // Information
                            Info(__FUNCTION__, "Selecting as numerator %s", objectN->GetName());

                            // Cast the object pointer into 2D histogram
                            TH2D * numerator2D = (TH2D*) processTH2(objectN);

                            // Numerator histogram
                            sprintf(name, "numerator_%s", denominator2D->GetName());
                            TH1D * numerator1D = numerator2D->ProjectionX(name, -1, -1, "e");

                            // MCTruth efficiencies histogram
                            sprintf(name, "mctruth_%s_%s", denominator2D->GetName(), numerator2D->GetName());
                            Info(__FUNCTION__, "Calculating efficiency %s", name);
                            TH1D * mctruth = (TH1D*) numerator1D->Clone();
                            mctruth->SetName(name);
                            mctruth->Divide(numerator1D, denominator1D, 1., 1., "e");
                            efficiencyHistogramSetup(mctruth);

                            // Save the efficiency histogram
                            output->cd("mctruth");
                            mctruth->Write();
                        }
                };
            }
    };

    if (!status) Error(__FUNCTION__, "Non matching histograms were found");

    return status;
}


bool PtrelTemplateMaker::makeTemplates(
    TFile * input,
    TFile * output
) const
{
    // Check function form exist
    for (std::size_t i = 0; i < (std::size_t)Flavor::Dimension; ++i)
        if (functions_[i].GetExpFormula().IsNull())
        {
            Error(__FUNCTION__, "Some or all function forms are not set");
            return false;
        }

    char name[256];

    // Information
    Info(__FUNCTION__, "Starting making templates");

    // Creating sub directory for templates
    output->cd();
    if (!output->FindKey("functions")) output->mkdir("functions");
    if (!output->FindKey("templates")) output->mkdir("templates");

    // Move to the directory with 2d histograms
    input->cd(PtrelTemplateMaker::directory);

    // Loop over all keys in this directory
    TIter nextkey( gDirectory->GetListOfKeys() );

    TKey * key;

    while (( key = (TKey*)nextkey() ))
    {
        TObject * object = key->ReadObj();
        if ( object->IsA()->InheritsFrom( "TH2" ) )
        {
            // Look over the flavor histograms
            for (Int_t i = 0; i < Flavor::Dimension; ++i)
                if ( containsIdentifier(object->GetName(),Flavor::Name[i]))
                {
                    // Cast the object pointer into 2D histogram
                    TH2D * histogram2D = (TH2D*) processTH2(object);

                    // Get the number of bins
                    Int_t nbins = histogram2D->GetNbinsX();

                    Double_t * lastParameters = 0;

                    for (Int_t j = 1; j <= nbins; ++j)
                    {
                        // Set histogram name by the formula = template_x_bin
                        sprintf(name, "template_%s_%d", histogram2D->GetName(), j);
                        // Project histogram
                        TH1D * histogram1D = histogram2D->ProjectionY(name, j, j, "e");

                        // Setup histogram
                        ptrelHistogramSetup(histogram1D);

                        // Clone a new function
                        TF1 * function = (TF1*) functions_[i].Clone();

                        // Set function name and title by the formula = function_x_bin
                        sprintf(name, "function_%s_%d", histogram2D->GetName(), j);
                        function->SetName(name);
                        function->SetTitle(name);

                        if (j != 1) function->SetParameters(lastParameters);

                        // Fit the histogram
                        Int_t fitStatus = histogram1D->Fit(function, "V", "", function->GetXmin(), function->GetXmax());

                        // Check the status of the fitting
                        if ( fitStatus )
                            Warning(__FUNCTION__, "Fitting problem returning status %d", fitStatus);
                        else
                            Info(__FUNCTION__, "Fitting %s chi2/ndf = (%f/%d)", histogram1D->GetName(), function->GetChisquare(), function->GetNDF());

                        // Get the parameter of the last optiomization
                        lastParameters = function->GetParameters();

                        // Saving the function and histograms
                        output->cd("functions");
                        function->Write();
                        output->cd("templates");
                        histogram1D->Write();
                    }
                }
        }
    };

    // Save rebinning information
    char factor[256];
    for (Int_t i = 0; i < Dependency::Dimension; ++i)
        if (rebin_[i] != 1)
        {
            output->cd();
            if (!output->FindKey("parameters")) output->mkdir("parameters");
            output->cd("parameters");
            Info(__FUNCTION__, "Saving rebinning factor %d for %s", rebin_[i], Dependency::Name[i]);
            sprintf(factor, "%d", rebin_[i]);
            TObjString(factor).Write(Dependency::Name[i]);
        }

    return true;
}
