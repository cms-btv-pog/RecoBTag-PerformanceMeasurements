#include "HistogramManager.h"
#include "system.h"

#include <TROOT.h>
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStyle.h"

ClassImp(HistogramManager)

HistogramManager::HistogramManager() 
{
 HistogramManager_debug = kFALSE ;
 
 postscriptfile = 0;
}

HistogramManager::~HistogramManager()
{
 if(postscriptfile) postscriptfile->Close();
}

////////////////////////////////////////////////////////////////////

void HistogramManager::WriteHistogram(TH1* anhisto, TString option)
{
 if(HistogramManager_debug) 
  this->Info("WriteHistogram(TH1* anhisto, Option_t)",
 	     "processing");
	     
 if(!anhisto) return ;
 
 anhisto->Write() ;
}

////////////////////////////////////////////////////
// macro to get an array of histo from a rootfile //
////////////////////////////////////////////////////

void HistogramManager::GetHistogramsInFile(TObjArray & table, TFile * afile)
{
 if(HistogramManager_debug) 
  this->Info("Get_Histos_in_File(TFile *afile, TObjArray & table)",
 	     "processing");
     
 if(!afile) { 
  this->Error("GetHistogramsInFile","no input file");
  return;
 } 
 
 table.Clear() ;
 
 Int_t goon = 0;
 TIter next(afile->GetListOfKeys());      
 TKey *key  = 0 ;

    while (goon < 3) {
       goon = 1;
       key = (TKey*)next();             
       if (!key) break;
       
       TObject *obj = key->ReadObj();   
       
       if (obj->InheritsFrom("TH1")) {  
//         if(HistogramManager_debug) this->"Found histo : " 
// 	                                << obj->GetName() << endl;
//         TH1 ahisto ;
//	 ahisto.Copy(*(TH1*)obj) ;
//	 table.Add((TH1*)ahisto.Clone());
	table.Add(obj);
       }
       else {
        this->Warning("Get_Histos_in_File(TFile *afile, TObjArray & table)",
		      "No histos were found in file"); 
       }
    goon++ ;
    }
}

///////////////////////////////////////////////////////////////////////////////

void HistogramManager::WriteAllHistogramsInFile(TString filename , 
						TString option)
{
 if(HistogramManager_debug) 
  this->Info("WriteAllHistogramsInFile(TString filename , TString option)",
 	     "processing");


 TObjArray histolist ;
 TString realfilename = filename ;
 TIter next(gDirectory->GetList()) ;
 TObject *obj = 0 ;
 
 Int_t goon = 0 ;

 while (goon < 2) {
       goon = 0;
       obj = (TObject*)next();  //get next obj in memory
       if (!obj) break;
       if (obj->InheritsFrom("TH1")) {  
        if(HistogramManager_debug) 
	 printf("Found histo : %s \n",obj->GetName());
	 histolist.Add(obj);
       }
    goon++;
    }   

 // Check if file already exists //

// if(FileExists(filename) && option.EndsWith("recreate")) {
// 
//  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
//  printf("!! Warning, file %s already exists while writing histos !!\n",filename.Data());
//  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
  
//  realfilename.Prepend("new_");

//  option = "recreate" ;
// }

 TFile file(realfilename,option);

 printf("==> Creating new file : %s to store histograms\n",file.GetName());
 
 for(Int_t i = 0; i < histolist.GetEntries(); i++)
  {
   TH1* histo = (TH1*)histolist[i];
   histo->Write();
  }
  file.Close(); 
}

/////////////////////////////////////////////////////////////////

void HistogramManager::SaveHistogram(TH1* hist1, Char_t* title , Char_t* option , 
	                             Bool_t logy , Bool_t logx , Int_t color ,
	                             Int_t stati , TString format)
{
 if(HistogramManager_debug) 
  this->Info("SaveHistogram()","processing");
	     
  if(format == "ps") {
    if(!postscriptfile) postscriptfile = new TPostScript("output.ps",111);
//    else postscriptfile->NewPage();
   }
  
 TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
 if (c1) c1->Clear();
 else c1 = new TCanvas("c1");
  
//  TPad* pad1 = new TPad("pad1","pad",0.03,0.03,0.98,0.98);
//  pad1->Draw();
  
  if (logy) c1->SetLogy(); 
  if (logx) c1->SetLogx();

//  pad1->cd();

  if (option == "COLZ") c1->Draw("SURF");
   
  gStyle->SetOptStat(stati);
  gStyle->SetPaperSize(20,28);
  
  hist1->SetTitle(title);

  hist1->SetFillColor(color); 

  hist1->GetXaxis()->SetTitleSize(0.05);
  hist1->Draw(option); 

// l'histo est memorise
  c1->Update();  
}

/////////////////////////////////////////////////////////////////////

void HistogramManager::SetHistogramStyle(TString style)
{
 if(style == "seb") SetSebStyle();
 //else if(style == "daniel") SetDanielStyle();
 else printf("Histogram style : %s does not exists\n",style.Data());
}

/////////////////////////////////////////////////////////////////////////

void HistogramManager::SetSebStyle(float W, float H, float titleW, float titleH)
{
// Graphics style parameters to avoid grey background on figures
   gStyle->SetCanvasColor(10);
   gStyle->SetStatColor(29);
   gStyle->SetTitleColor(42);
   gStyle->SetPadColor(10);
   gStyle->SetPaperSize(20,24);
   gStyle->SetStatFont(0);
   gStyle->SetOptDate(0);

   //title size

   gStyle->SetOptStat(111111);
   gStyle->SetOptFit(0);
   gStyle->SetStatW(W); 
   gStyle->SetStatH(H); 
   gStyle->SetTitleFont(1);
   if (titleH>0) gStyle->SetTitleH(titleH);
   if (titleW>0) gStyle->SetTitleW(titleW); 
   gStyle->SetFillColor(2); 
   gStyle->SetMarkerColor(2); 
 
   gStyle->SetMarkerStyle(8); 
   gStyle->SetMarkerSize(0.3);
   
   gStyle->SetLabelSize(0.05);
   gStyle->SetLabelSize(0.05);

   gStyle->SetTitleOffset(1.3);
   gStyle->SetTitleSize(0.05);
}
