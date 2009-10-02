#include <iostream>

#include "TString.h"
#include "TSystem.h"
#include "TFile.h"

#include <fstream>

TFile   *fnew;
TList   *flist;
TFile   *afile, *file1;

TH1     *h1, *h2;
TObject *obj;
TKey    *key;

Bool_t debug = false ;

void AddRecursive(TDirectory *root,TDirectory* node);

void AddHist(TString list_of_files, TString outputfilename)
{
// Create a new output file //
fnew = new TFile(outputfilename.Data(),"new") ;

//create a support list for the input files
flist = new TList();

//open all input files and insert them in the list of files
Int_t nfiles = 0;

ifstream list(list_of_files.Data(),ios::in);
if(!list) {
cout << "!! WARNING : no list : " << list_of_files << " found in current directory" << endl; return ;
}
 
  while (!list.eof()) {
  
  TString filename = "" ;
  
  list >> filename ;
  if(filename == "") break ;
  
    afile = new TFile(filename.Data());
    cout << "Adding file : " << filename <<  " to the list of files" << endl;
    flist->Add(afile);
    nfiles++;
  }

//Get a pointer to the first file
afile = file1 = (TFile*)flist->First();

 AddRecursive(fnew,file1);

if(debug) fnew->ls();
  fnew->Write();
  fnew->Close();
  delete fnew;
  flist->Delete();
  delete flist;
}

//////////////////////////////////////////////////////////////////////////////

void AddRecursive(TDirectory *root,TDirectory* node) {

  static TDirectory *dact;

  TDirectory *dirsav;

//We create an iterator to loop on all objects(keys) of first file
TIter nextkey(node->GetListOfKeys());

  while (key = (TKey*)nextkey()) {
    node->cd();
    obj = key->ReadObj();
 
    if(obj->InheritsFrom("TH1")) { 
    
    if(debug) cout << "Getting histo : " << obj->GetName() << endl;
    
      h1 = (TH1*)obj;
      afile = (TFile*)flist->After(file1);
      
      while (afile) { //loop on all files starting at second file
        char* base=strstr(root->GetPath(),":"); base+=2;

        dirsav = gDirectory;
        afile->cd(base);
        h2 = (TH1*)gDirectory->Get(h1->GetName());
	if(debug) cout << "Adding histo : " << h2->GetName() << " from file : " <<
	gDirectory->GetName() << endl;
        
	dirsav->cd();
        if (h2) { 
	
         h1->Add(h2);
          delete h2;
        }
        afile = (TFile*)flist->After(afile);
      }
    } 

// write node object, modified or not into fnew
if (obj) {
      root->cd();
      obj->Write(key->GetName());
      delete obj;
      obj=NULL;
    }
  }
  root->cd();
}




