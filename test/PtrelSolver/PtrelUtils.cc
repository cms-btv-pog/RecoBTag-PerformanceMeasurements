
#include <cmath>
#include <iostream>

#include "PtrelUtils.h"
#include "PtrelSolverDependencies.h"

void ptrelHistogramSetup(TH1 * histogram)
{
    // Set histogram name
    histogram->SetTitle(histogram->GetName());
    // Set histogram x title
    histogram->GetXaxis()->SetTitle("p_{T}^{rel} [GeV/c]");
    // Set histogram x title
    histogram->GetYaxis()->SetTitle("counts");
}

void efficiencyHistogramSetup(TH1 * histogram)
{
    // Set histogram name
    histogram->SetTitle(histogram->GetName());

    for (Int_t i = 0; i < Dependency::Dimension; ++i)
        if (containsIdentifier(histogram->GetName(),Dependency::Name[i]))
        {
            // Set histogram x title
            histogram->GetXaxis()->SetTitle(Dependency::Label[i]);
            break;
        }

    for (Int_t i = 0; i < Flavor::Dimension; ++i)
        if (containsIdentifier(histogram->GetName(), Flavor::Name[i]))
        {
            char name[256];
            sprintf(name, "%s efficiency", Flavor::Label[i]);
            histogram->GetYaxis()->SetTitle(name);
            break;
        }
}

// Look for an identifier in histogram name
bool containsIdentifier(const char* objName_, const char* id_)
{
    bool ret=false;
    TString objName(objName_);
    TString id(id_);
    TString p1;
    p1+="_";
    p1+=id;
    p1+="_";
    TString p2;
    p2+="_";
    p2+=id;
    p2+="$";
    TPRegexp reg1(p1);
    TPRegexp reg2(p2);
    if (objName.Contains(reg1) || objName.Contains(reg2)) ret=true;
    return ret;
}


void ErrorPropagator::transformation(const TMatrixD & matrix)
{
    // Save the matrix
    matrix_.ResizeTo(matrix);
    matrix_ = matrix;
    matrix_.Print();

    // Define a place holder for the Hadamard product
    hadamard_.ResizeTo(matrix);
    inverseHadamard_.ResizeTo(matrix);

    // Compute the Hadamard product
    for (Int_t i = 0; i < matrix.GetNrows(); ++i)
        for (Int_t j = 0; j < matrix.GetNcols(); ++j)
        {
            hadamard_(i,j) = matrix(i,j) * matrix(i,j);
            inverseHadamard_(i,j) = hadamard_(i,j);
        }

    // Compute the inverse of Hadamard matrix
    inverseHadamard_.Invert();
    inverseHadamard_.Print();
}

TVectorD const & ErrorPropagator::correctedErrors(TVectorD const & errors)
{
    correctedErrors_.ResizeTo(errors.GetNrows());

    if (removeCorrelations_)
    {
        correctedErrors_ = errors;
        return correctedErrors_;
    }

    TMatrixD holder(errors.GetNrows(), 1);
    for (Int_t k = 0; k < errors.GetNrows(); ++k)
        holder(k,0) = errors(k);

    TMatrixD noCorrelatedErrors(inverseHadamard_, TMatrixD::kMult, holder);

    bool correctedErrorFlag = false;

    for (Int_t k = 0; k < noCorrelatedErrors.GetNrows(); ++k)
        if (noCorrelatedErrors(k, 0) < 0.)
        {
            noCorrelatedErrors(k, 0) = 0.;
            correctedErrorFlag = true;
        }

    if (correctedErrorFlag)
    {
        std::cout << "Warning in <ErrorPropagator::correctedErrors()>: ";
        std::cout << "some errors are not consistent with the correlations (they will be corrected)." << std::endl;
    }

    TMatrixD correctedErrors(hadamard_, TMatrixD::kMult, noCorrelatedErrors);

    for (Int_t k = 0; k < correctedErrors.GetNrows(); ++k)
        correctedErrors_(k) = correctedErrors(k,0);

    return correctedErrors_;
}


Double_t ErrorPropagator::operator()(TVectorD const & derivatives, TVectorD const & errors)
{
    correctedErrors(errors);

    Double_t result = 0;

    for (Int_t k = 0; k < correctedErrors_.GetNrows(); ++k)
    {
        result += derivatives(k) * derivatives(k) * correctedErrors_(k) * correctedErrors_(k);

        if (!removeCorrelations_)
        {
            for (Int_t i = 0; i < derivatives.GetNrows(); ++i)
                for (Int_t j=i+1; j < derivatives.GetNrows(); ++j)
                {
                    Double_t element = 0;

                    for (Int_t n = 0; n < errors.GetNrows(); ++n)
                        element += matrix_(i,n)*matrix_(j,n)*inverseHadamard_(n,k);

                    result += 2 * element * derivatives(i) * derivatives(j) * correctedErrors_(k) * correctedErrors_(k);
                }
        }
    }

    return sqrt(result);
}

