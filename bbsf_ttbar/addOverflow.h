#include <TH1D.h>
#include <iostream>

void addOverflow(TH1D * h)
{
   const int n = h->GetNbinsX();
   std::cout << "entries in overflow: " << h->GetBinContent(n+1) << " from hist " << TString(h->GetName()) + "..." << std::endl;
   h->AddBinContent(n, h->GetBinContent(n+1));
   h->SetBinContent(n+1, 0);
}

