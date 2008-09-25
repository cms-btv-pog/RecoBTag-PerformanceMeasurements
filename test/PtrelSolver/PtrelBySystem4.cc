
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

#include "PtrelUtils.h"

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

        // Collection of flavor content n in sample
        ValueVector pValues;
        CovarianceVector pCovariances;

        // Name for the p sample
        sprintf(name, "%s/p_%s", directory, Dependency::Name[i]);
        sprintf(nsample, "p_%s", Dependency::Name[i]);
        // Measure the flavor content in n sample
        CallSafely( measure(input, output, name, pValues, pCovariances) )
  
        // Evauluate the coefficients between b-efficiencies from p and n samples
        sprintf(name, "mctruth_n_%s_b_ntag_%s_b_[A-Z]*$", Dependency::Name[i], Dependency::Name[i]);
        TPRegexp defficiency(name);
        sprintf(name, "mctruth_p_%s_b_ptag_%s_b_[A-Z]*$", Dependency::Name[i], Dependency::Name[i]);
        TPRegexp nefficiency(name);        
        CallSafely( makeCoefficients(templates_, output, defficiency, nefficiency) )

        // Evauluate the coefficients between b-efficiencies from p and n samples
        sprintf(name, "mctruth_n_%s_cl_ntag_%s_cl_[A-Z]*$", Dependency::Name[i], Dependency::Name[i]);
        TPRegexp defficiency2(name);
        sprintf(name, "mctruth_p_%s_cl_ptag_%s_cl_[A-Z]*$", Dependency::Name[i], Dependency::Name[i]);
        TPRegexp nefficiency2(name);        
        CallSafely( makeCoefficients(templates_, output, defficiency2, nefficiency2) )
        
        // HACK to force reading the objets from file
        // Close output
        input->Close();
        // Reopen
        GetSafely(input, TFile::Open(inputfile, "READ"))        
        // ENDHACK

        // Get a 2D histogram to produce a place holder for the results later
        sprintf(name, "%s/n_%s", directory, Dependency::Name[i]);
        CreateSafely(TH2D, histogram2D, readTH2(input, name))

        /*for (std::size_t j = 0; j < ptagHistograms.size(); ++j)
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
        } */
    }

    input->Close();
    output->Close();
}


bool PtrelBySystem4::makeCoefficients (
    TFile * input,
    TFile * output,
    TPRegexp & patternD,
    TPRegexp & patternN
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
    input->cd("/mctruth/");

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


/*bool PtrelBySystem4::compute(
    TH1 * histogram,
    ValueVector const & factors,
    TVectorD const & ntagCounters,
    TVectorD const & ptagCounters,
    ValueVector const & nValues,
    CovarianceVector const & nCovariance,
    ValueVector const & pValues,
    CovarianceVector const & pCovariance
)
{
    // Loop over different bins
    for (std::size_t i = 0; i < nValues.size(); ++i)
    {
    	Double_t n_tag = ntagCounters[i];
        Double_t n_b = nValues[i](Flavor::b);
        Double_t n_cl = nValues[i](Flavor::cl);
    	Double_t p_tag = ptagCounters[i];
        Double_t p_b = ptagValues[i](Flavor::b); 
        Double_t p_cl = ptagValues[i](Flavor::cl); 
        
        // Double_t n_cl_error = nCovariance[i](Flavor::cl, Flavor::cl);
        // Double_t p_b_error = ptagCovariance[i](Flavor::b, Flavor::b);

        Double_t alfa = factors[i][Flavor::cl];
        Double_t beta = factors[i][Flavor::b];

        histogram->SetBinContent(i+1, (alfa*n_tag*p_cl - p_tag*n_cl)/(alfa*n_b*p_cl - beta*n_cl*p_b));
    }
    return true;
}
*/
