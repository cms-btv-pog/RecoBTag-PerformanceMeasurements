
#include <iostream>

#include "BasePlotter.h"


ClassImp(BasePlotter)


Long64_t BasePlotter::LoadTree(Long64_t entry)
{
  if (!fChain_) return -5;

  Long64_t centry = fChain_->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain_->IsA() != TChain::Class()) return centry;

  TChain *chain = (TChain*)fChain_;
  if (chain->GetTreeNumber() != fCurrent_)
    fCurrent_ = chain->GetTreeNumber();

  return centry;
}


void BasePlotter::Loop(Long64_t max)
{
  if (fChain_ == 0) return;

  // Reset the current counter
  fCurrent_ = -1;

  // Set the branch address
  fChain_->SetBranchAddress("s8.", &event_);

  // Book histograms
  Book();
  
  Long64_t nentries = fChain_->GetEntriesFast();
  std::cout << " Total entries = " << fChain_->GetEntries() << std::endl;
  
  // Loop over the events in the chain
  Long64_t nbytes = 0;
  
  for (Long64_t jentry=0; jentry<nentries; jentry++)
  {
  	if (max != 0 && jentry >= max) break;  
  	
    Long64_t ientry = LoadTree(jentry);
    
    if (ientry < 0) break;
  	
    nbytes += GetEntry(jentry);   
  	
    if((jentry+1)%10000 == 0) 
      std::cout << "### processing entry : " << (jentry+1) << " (" << ((double)nbytes/1000000) << " Mb)" << std::endl;

    Fill(event_);
  }  
}

