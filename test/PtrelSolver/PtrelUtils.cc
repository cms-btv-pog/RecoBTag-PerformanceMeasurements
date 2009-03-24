
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
bool containsIdentifier(const char* objName_, const char* id_) {
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
    if(objName.Contains(reg1) || objName.Contains(reg2)) ret=true;
    return ret;
}

