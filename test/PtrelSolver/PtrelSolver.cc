#include "PtrelSolver.h"

#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TMinuit.h"
#include "TObjString.h"

#include "PtrelUtils.h"


ClassImp(PtrelSolver)


char const * PtrelSolver::directory = "/Histograms/muon_in_jet";


void PtrelSolver::templates(const char * filename)
{
    if (templates_) templates_->Close();

    GetSafely(templates_, TFile::Open(filename, "READ"))

    char name[256];

    for (Int_t i = 0; i < Dependency::Dimension; ++i)
    {
        sprintf(name, "/parameters/%s", Dependency::Name[i]);
        CreateSafely(TObjString, value, templates_->Get(name))
        rebin_[i] = value->String().Atoi();
        Info(__FUNCTION__, "Loading rebinning factor %d for %s", rebin_[i], Dependency::Name[i]);
    }
}


bool PtrelSolver::measure(
    TFile * input,
    TFile * output,
    const char * name,
    ValueVector & vVector,
    CovarianceVector & cVector
)
{
    CreateSafelyZero(TH2D, histogram, readTH2(input, name))
    CallSafelyZero( measure(output, histogram, vVector, cVector) )
    return true;
}


bool PtrelSolver::measure(
    TFile * input,
    TFile * output,
    TPRegexp pattern,
    StringVector & hVector,
    ValueMatrix & vMatrix,
    CovarianceMatrix & cMatrix
)
{
    // Clean the containers
    hVector.clear();
    vMatrix.clear();
    cMatrix.clear();

    // Move to the directory with 2d histograms
    input->cd(directory);

    // Loop over all keys in this directory
    TIter nextkey( gDirectory->GetListOfKeys() );

    TKey * key;

    //char name[256];

    while ( key = (TKey*)nextkey() )
    {
        // Measure those 2D histogram with name fit the pattern
        TObject * object = key->ReadObj();
        if ( object->IsA()->InheritsFrom( "TH2" ) )
            if ( TString(object->GetName() ).Contains(pattern) )
            {
                // Temporal container for efficiencies
                TH2 * histogram;
                ValueVector values;
                CovarianceVector covariances;

                // Process 2D histogram associated to the keyword
                GetSafelyZero(histogram, processTH2(object))

                // Measuring efficiencies
                CallSafelyZero( measure(output, histogram, values, covariances) )

                // Appending the results
                hVector.push_back( TString(histogram->GetName()) );
                cMatrix.push_back(covariances);
                vMatrix.push_back(values);
            }
    };

    return true;
}


bool PtrelSolver::measure(
    TFile * output,
    TH2 * histogram2D,
    ValueVector & vVector,
    CovarianceVector & cVector
)
{
    char name[256];

    // Creating sub directory for templates
    if (!output->FindKey("fits")) output->mkdir("fits");

    // Clearing the containers
    vVector.clear();
    cVector.clear();

    // Number of bins
    Int_t nbins = histogram2D->GetNbinsX();

    // Loop over the bins
    for (Int_t i = 1; i <= nbins; ++i)
    {
        // Information
        Info(__FUNCTION__, "Measuring flavor content in bin %d", i);

        // Get the proper combined function
        CreateSafelyZero(TF1, function, combinedFunction(histogram2D->GetName(), i))

        // Set histogram name by the formula = fit_x_bin
        sprintf(name, "fit_%s_%d", histogram2D->GetName(), i);
        // Project histogram
        TH1D * histogram1D = histogram2D->ProjectionY(name, i, i, "e");
        // Basic setup
        ptrelHistogramSetup(histogram1D);

        // Temporal container with the results
        TVectorD values;
        TMatrixD covariance;

        // Fitting the histogram
        CallSafelyZero( fit(histogram1D, function, values, covariance) )

        // Collecting the results
        vVector.push_back(values);
        cVector.push_back(covariance);

        // Saving the fitting.
        output->cd();
        output->cd("fits");
        histogram1D->Write();
    }
    return true;
}


TH2 * PtrelSolver::readTH2(TFile * file, const char * name) const
{
    CreateSafelyZero(TObject, object, file->Get(name))
    CreateSafelyZero(TH2, histogram, processTH2(object))
    return histogram;
}


TH2 * PtrelSolver::processTH2(TObject * object) const
{
    // Cast the object pointer into 2D histogram
    TH2 * histogram = (TH2*) object;

    // Rebinning in ptrel
    if ( rebin_[Dependency::ptrel] != 1 )
    {
        Info(__FUNCTION__, "Rebinning ptrel a factor %d in %s", rebin_[Dependency::ptrel], histogram->GetName());
        histogram->RebinY(rebin_[Dependency::ptrel]);
    }

    // Check the histogram dependency and rebin
    for (Int_t i = 1; i < Dependency::Dimension; ++i)
        if ( TString(histogram->GetName()).Contains(Dependency::Name[i]) && rebin_[i] != 1)
        {
            Info(__FUNCTION__, "Rebinning %s a factor %d in %s", Dependency::Name[i], rebin_[i], histogram->GetName());
            histogram->RebinX(rebin_[i]);
        }

    return histogram;
}


bool PtrelSolver::fit(TH1 * histogram, TF1 * function, TVectorD & values, TMatrixD & covariance)
{
    // Check for initial values (otherwise using default)
    if (values.GetNoElements() != Flavor::Dimension)
    {
        Warning(__FUNCTION__, "Inconsistent dimension for initial values using default");
        values.ResizeTo(Flavor::Dimension);
        values.Zero();
        values(1) = 1.;
    }

    // Resizing error matrix
    covariance.ResizeTo(Flavor::Dimension, Flavor::Dimension);

    // Total number of event including correction because of weights
    Double_t numberEvents = histogram->Integral();

    // Calculation of normalization parameters
    for (Int_t i = 0; i < Flavor::Dimension; ++i)
    {
        function->SetParameter(i, numberEvents * values(i));
        Info(__FUNCTION__, "Fitting with initial fraction %f for %s", values(i), Flavor::Name[i]);
    }

    // Fitting combined function against data
    histogram->Fit(function, "Q", "", function->GetXmin(), function->GetXmax());

    // Information
    Info(__FUNCTION__, "Fitting %s chi2/ndf = (%f/%d)", histogram->GetName(), function->GetChisquare(), function->GetNDF());

    // Temporal containers for the results
    TVectorD values_(Flavor::Dimension);
    TMatrixD covariance_(Flavor::Dimension, Flavor::Dimension);

    // Getting the results for the fitting
    function->GetParameters(values_.GetMatrixArray());
    // Get the covariance matrix of parameter errors
    gMinuit->mnemat(covariance_.GetMatrixArray(), Flavor::Dimension);

    // Copy semantics
    values = values_;
    covariance = covariance_;

    // Scale to the number of events
    Double_t scale = numberEvents/values.Sum();
    values *= scale;
    covariance *= (scale*scale);

    return true;
}


TF1 * PtrelSolver::combinedFunction(char const * keyword, Int_t bin)
{
    char name[256];
    TString form;
    Double_t xmin = 0;
    Double_t xmax = 0;

    // Look over the flavor functions and creating
    for (Int_t i = 0; i < Flavor::Dimension; ++i)
    {
        // HACK to solve the problem of naming conventio
        TPRegexp p("[A-Z]{3,}$");
        if (TString(keyword).Contains(p))
        {
            TPRegexp t("[A-Z]{3,}");
            Int_t inx = TString(keyword).Index(t);
            std::string tmp(keyword);
            std::string tag(tmp.substr(inx, tmp.size()));
            std::string core(tmp.substr(0, inx-1));
            // Read the function from the file
            sprintf(name, "/functions/function_%s_%s_%s_%d", core.c_str(), Flavor::Name[i], tag.c_str(), bin);
        }  // END HACK
        else
            sprintf(name, "/functions/function_%s_%s_%d", keyword, Flavor::Name[i], bin);

        Info(__FUNCTION__, "Loading %s", name);
        CreateSafelyZero(TF1, function, templates_->Get(name))

        // Collect minimal and maximal values
        xmin = function->GetXmin();
        xmax = function->GetXmax();

        // Calculate the integral
        Double_t scale = 1./function->Integral(xmin, xmax);

        // Create the funcion form for combined function
        if (i != 0) form += '+';
        form += '[';
        form += i;
        form += "]*";
        form += scale;
        form += "*(";
        form += function->GetExpFormula("p");
        form += ')';
    }

    // Information
    sprintf(name, "combined_%s_%d", keyword, bin);
    Info(__FUNCTION__, "Creating combined function %s", name);

    // Creating function
    return new TF1(name, form.Data(), xmin, xmax);
}
