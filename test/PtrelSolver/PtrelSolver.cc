
#include <math.h>

#include "PtrelSolver.h"

#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TMinuit.h"
#include "TObjString.h"

#include "RecoBTag/PerformanceMeasurements/interface/CFractionFitter.h"

#include "PtrelUtils.h"


ClassImp(PtrelSolver)


char const * PtrelSolver::directory = "/muon_in_jet";


void PtrelSolver::templates(const char * filename)
{
    if (templates_) templates_->Close();

    GetSafely(templates_, TFile::Open(filename, "READ"))

    char name[256];

    for (Int_t i = 0; i < Dependency::Dimension; ++i)
    {
        sprintf(name, "/parameters/%s", Dependency::Name[i]);
        TObjString * value = (TObjString *) templates_->Get(name);
        if (value)
        {
            rebin_[i] = value->String().Atoi();
            Info(__FUNCTION__, "Loading rebinning factor %d for %s", rebin_[i], Dependency::Name[i]);
        }
    }
}


bool PtrelSolver::measure(
    TFile * input,
    TFile * output,
    const char * name,
    ValueVector & vVector,
    ValueVector & eVector
)
{
    CreateSafelyZero(TH2D, histogram, readTH2(input, name))
    CallSafelyZero(measure(output, histogram, vVector, eVector))
    return true;
}


bool PtrelSolver::measure(
    TFile * input,
    TFile * output,
    TPRegexp pattern,
    StringVector & hVector,
    ValueMatrix & vMatrix,
    ValueMatrix & eMatrix
)
{
    // Return status
    bool status = false;

    // Clean the containers
    hVector.clear();
    vMatrix.clear();
    eMatrix.clear();

    // Move to the directory with 2d histograms
    input->cd(directory);

    // Loop over all keys in this directory
    TIter nextkey( gDirectory->GetListOfKeys() );

    TKey * key;

    //char name[256];

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
                TH2 * histogram;
                ValueVector values, errors;

                // Process 2D histogram associated to the keyword
                GetSafelyZero(histogram, processTH2(object))

                // Measuring efficiencies
                CallSafelyZero( measure(output, histogram, values, errors) )

                // Appending the results
                hVector.push_back(TString(histogram->GetName()));
                vMatrix.push_back(values);
                eMatrix.push_back(errors);
            }
    };

    if (!status) Error(__FUNCTION__, "Non matching histograms were found");

    return status;
}


bool PtrelSolver::measure(
    TFile * output,
    TH2 * histogram2D,
    ValueVector & vVector,
    ValueVector & eVector
)
{
    char name[256];

    // Creating sub directory for templates
    if (!output->FindKey("fits")) output->mkdir("fits");

    // Clearing the containers
    vVector.clear();
    eVector.clear();

    // Number of bins
    Int_t nbins = histogram2D->GetNbinsX();

    // Loop over the bins
    for (Int_t i = 1; i <= nbins; ++i)
    {
        // Information
        Info(__FUNCTION__, "Measuring flavor content using fit option %s in bin %d", Fit::Label[fittype_], i);

        // Get the proper set of histograms
        CallSafelyZero( combinedTemplates(histogram2D->GetName(), i) )

        // Temporal container with the results
        TVectorD values, errors;

        if (fittype_ == Fit::histograms)
        {
            // Set histogram name by the formula = data_x_bin (transient)
            sprintf(name, "data_%s_%d", histogram2D->GetName(), i);
            // Project histogram
            TH1D * histogram1D = histogram2D->ProjectionY(name, i, i, "e");

            // Fitting the histogram
            CallSafelyZero(fit(output, histogram1D, combinedHistograms_, values, errors))
        }
        else
        {
            // Set histogram name by the formula = fit_x_bin
            sprintf(name, "fit_%s_%d", histogram2D->GetName(), i);
            // Project histogram
            TH1D * histogram1D = histogram2D->ProjectionY(name, i, i, "e");
            // Basic setup
            ptrelHistogramSetup(histogram1D);

            // Fitting the histogram
            CallSafelyZero(fit(histogram1D, combinedFunctions_, values, errors))

            // Saving the fitting.
            output->cd();
            output->cd("fits");
            histogram1D->Write();
        }
        // Collecting the results
        vVector.push_back(values);
        eVector.push_back(errors);
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
        if ( containsIdentifier(histogram->GetName(),Dependency::Name[i]) && rebin_[i] != 1)
        {
            Info(__FUNCTION__, "Rebinning %s a factor %d in %s", Dependency::Name[i], rebin_[i], histogram->GetName());
            histogram->RebinX(rebin_[i]);
        }

    return histogram;
}


bool PtrelSolver::fit(TH1 * histogram, TF1 * function, TVectorD & values, TVectorD & errors)
{
    // Check the number of flavor is not lower than two, otherwise cancel the fitting.
    if ( fitFlavors_.size() < 2 )
    {
        Error(__FUNCTION__, "Number of flavor to fit is lower than two.");
        return false;
    }

    // Check for initial values (otherwise using default)
    if (values.GetNoElements() != (Int_t) fitFlavors_.size())
    {
        Warning(__FUNCTION__, "Inconsistent dimension for initial values using default");
        values.ResizeTo(fitFlavors_.size());
        values.Zero();
        for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
            values(i) = 1./fitFlavors_.size();
    }

    // Numerical integration of data histogram in the interval of the function
    Double_t integral = histogram->Integral(
                            histogram->FindBin(function->GetXmin()),
                            histogram->FindBin(function->GetXmax()),
                            "width"
                        );

    // Calculation of normalization parameters
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
    {
        function->SetParameter(i, integral * values(i));
        Info(__FUNCTION__, "Fitting with initial fraction %f for %s", values(i), Flavor::Name[fitFlavors_[i]]);
    }

    // Fitting combined function against data
    histogram->Fit(function, "Q", "", function->GetXmin(), function->GetXmax());

    // Information
    Info(__FUNCTION__, "Fitting %s chi2/ndf = (%f/%d)", histogram->GetName(), function->GetChisquare(), function->GetNDF());

    // Temporal containers for the results
    TVectorD values_(fitFlavors_.size());
    TMatrixD covariance_(fitFlavors_.size(), fitFlavors_.size());

    // Getting the results for the fitting
    function->GetParameters(values_.GetMatrixArray());
    // Get the covariance matrix of parameter errors
    gMinuit->mnemat(covariance_.GetMatrixArray(), fitFlavors_.size());

    // Copy semantics for the values
    values = values_;

    // Total number of event including correction because of weights
    Double_t numberEvents = histogram->Integral();

    // Sum of the x values
    Double_t sumx = values.Sum();

    // Derivative matrix for error computation
    TMatrixD derivative(fitFlavors_.size(), fitFlavors_.size());
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
        for (Int_t j = 0; j < (Int_t) fitFlavors_.size(); ++j)
        {
            if (i == j)
                derivative(j,i) = 1./sumx - values(i)/(sumx*sumx);
            else
                derivative(j,i) = - values(i)/(sumx*sumx);
        }

    // Error computation
    TMatrixD temporal (covariance_, TMatrixD::kMult, derivative);
    TMatrixD covariance (derivative.T(), TMatrixD::kMult, temporal);

    errors.ResizeTo(fitFlavors_.size());
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
        errors(i) = sqrt( numberEvents*numberEvents*covariance(i,i) + numberEvents*values(i)*values(i)/(sumx*sumx) );

    // Scaling to the number of events
    values *= numberEvents/sumx;

    return true;
}


bool PtrelSolver::fit(TFile * output, TH1 * histogram, TObjArray * templates, TVectorD & values, TVectorD & errors)
{
    // Check the number of flavor is not lower than two, otherwise cancel the fitting.
    if ( fitFlavors_.size() < 2 )
    {
        Error(__FUNCTION__, "Number of flavor to fit is lower than two.");
        return false;
    }

    // Check if data histogram is empty
    if (!histogram->GetEntries())
    {
        Warning(__FUNCTION__, "Empty data histogram !");
        return false;
    }

    // Creating a fraction fitter object
    CFractionFitter fit(histogram, templates);

    // Set the constrains for all fraction values.
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
        fit.Constrain(i, -1.0, 2.0);

    // Running the fitter
    Int_t status = fit.Fit();

    // Resizing error vector
    values.ResizeTo(fitFlavors_.size());
    errors.ResizeTo(fitFlavors_.size());

    // Warning if the fit do not converge
    if (status)
    {
        Warning(__FUNCTION__, "Fitting do not converge (status code %d)", status);
        values.Zero();
        errors.Zero();
        return true;
    }

    // Information
    Info(__FUNCTION__, "Fitting %s chi2/ndf = (%f/%d)", histogram->GetName(), fit.GetChisquare(), fit.GetNDF());

    // Resizing error vector
    values.ResizeTo(fitFlavors_.size());
    errors.ResizeTo(fitFlavors_.size());

    // Total number of event including correction because of weights
    Double_t numberEvents = histogram->Integral();

    // Getting the results and scaling up to the number of events
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
    {
        Double_t value, error;
        fit.GetResult(i, value, error);
        values(i) = numberEvents * value;
        errors(i) = sqrt( numberEvents*value*value + numberEvents*numberEvents*error*error );
    }

    // Saving the fitting.
    output->cd();
    output->cd("fits");
    TH1 * result = fit.GetPlot();
    ptrelHistogramSetup(result);
    result->Write();

    return true;
}


bool PtrelSolver::combinedTemplates(char const * keyword, Int_t bin)
{
    char functionName[256];
    char templateName[256];

    TString form;
    Double_t xmin = 0;
    Double_t xmax = 0;

    combinedHistograms_ = new TObjArray(fitFlavors_.size());

    // Look over the flavor functions and creating
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
    {
        // Reading tagged template funtions
        TPRegexp p("[A-Z]{3,}$");
        if (TString(keyword).Contains(p))
        {
            TPRegexp t("[A-Z]{3,}");
            Int_t inx = TString(keyword).Index(t);
            std::string tmp(keyword);
            std::string tag(tmp.substr(inx, tmp.size()));
            std::string core(tmp.substr(0, inx-1));

            // Function name to be got from the file
            if (!taggedLightTemplate_ && fitFlavors_[i] == Flavor::l)
            {
                TPRegexp p("tag");
                inx = TString(core).Index(p);
                std::string dependency(core.substr(inx+4,core.size()));
                sprintf(functionName, "/functions/function_n_%s_%s_%d", dependency.c_str(), Flavor::Name[fitFlavors_[i]], bin);
                sprintf(templateName, "/templates/template_n_%s_%s_%d", dependency.c_str(), Flavor::Name[fitFlavors_[i]], bin);
            }
            else
            {
                sprintf(functionName, "/functions/function_%s_%s_%s_%d", core.c_str(), Flavor::Name[fitFlavors_[i]], tag.c_str(), bin);
                sprintf(templateName, "/templates/template_%s_%s_%s_%d", core.c_str(), Flavor::Name[fitFlavors_[i]], tag.c_str(), bin);
            }
        }
        else
        {
            // Function name to be got from the file
            if (!taggedLightTemplate_ && fitFlavors_[i] == Flavor::l)
            {
                TPRegexp p("(n|p)_");
                Int_t inx = TString(keyword).Index(p);
                std::string tmp(keyword);
                std::string dependency(tmp.substr(inx+2, tmp.size()));
                sprintf(functionName, "/functions/function_n_%s_%s_%d", dependency.c_str(), Flavor::Name[fitFlavors_[i]], bin);
                sprintf(templateName, "/templates/template_n_%s_%s_%d", dependency.c_str(), Flavor::Name[fitFlavors_[i]], bin);
            }
            else
            {
                sprintf(functionName, "/functions/function_%s_%s_%d", keyword, Flavor::Name[fitFlavors_[i]], bin);
                sprintf(templateName, "/templates/template_%s_%s_%d", keyword, Flavor::Name[fitFlavors_[i]], bin);
            }
        }

        Info(__FUNCTION__, "Loading %s", functionName);
        CreateSafelyZero(TF1, function, templates_->Get(functionName))

        Info(__FUNCTION__, "Loading %s", templateName);
        CreateSafelyZero(TH1, histogram, templates_->Get(templateName))

        // Collect minimal and maximal values
        xmin = function->GetXmin();
        xmax = function->GetXmax();

        Double_t sum = 0;

        // Sum of the function over the bin centers
        for (Int_t j=1; j<= histogram->GetNbinsX(); ++j)
            sum += function->Eval(histogram->GetBinCenter(j));

        // Scale the function to sum = 1
        Double_t scale = 1./sum;

        // Create the funcion form for combined function
        if (i != 0) form += '+';
        form += '[';
        form += i;
        form += "]*";
        form += scale;
        form += "*(";
        form += function->GetExpFormula("p");
        form += ')';

        // Add the histograms
        if (histogram->GetEntries() == 0)
        {
            Error(__FUNCTION__, "Empty histogram %s", templateName);
            return false;
        }
        combinedHistograms_->Add(histogram);
    }

    // Information
    sprintf(functionName, "combined_%s_%d", keyword, bin);
    Info(__FUNCTION__, "Creating combined function %s defined as %s", functionName, form.Data());

    // Creating function
    combinedFunctions_ = new TF1(functionName, form.Data(), xmin, xmax);

    return true;
}

