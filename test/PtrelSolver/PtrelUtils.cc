
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
        if (TString(histogram->GetName()).Contains(Dependency::Name[i]))
        {
            // Set histogram x title
            histogram->GetXaxis()->SetTitle(Dependency::Label[i]);
            break;
        }

    for (Int_t i = 0; i < Flavor::Dimension; ++i)
        if (TString(histogram->GetName()).Contains(Flavor::Name[i]))
        {
            char name[256];
            sprintf(name, "%s efficiency", Flavor::Label[i]);
            histogram->GetYaxis()->SetTitle(name);
            break;
        }
}

