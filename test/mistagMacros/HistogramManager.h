#ifndef HistogramManager_h
#define HistogramMangaer_h

#include "TString.h"
#include "TFile.h"
#include "TPostScript.h"

#include "TH2F.h"
#include "TH2D.h"
#include "TH2C.h"
#include "TH2S.h"

#include "TH3F.h"
#include "TH3D.h"
#include "TH3C.h"
#include "TH3S.h"

#include "TProfile.h"

#include "TObjArray.h"

class HistogramManager : public TObject 
{
 private :

 Bool_t HistogramManager_debug ;
 
  // Rajouter un TFile* //
  // une methode qui initialise ce TFile (on peut en ecrire plusiseurs)
  // 
  
 TPostScript* postscriptfile ;

 public:

 // Accessible method //

 HistogramManager();

 virtual ~HistogramManager();

 void SetDebug(Bool_t adebug) { HistogramManager_debug = adebug; }
 
 void SetPostScriptFile(TPostScript & apsfile) { postscriptfile = &apsfile;}

 void GetHistogramsInFile(TObjArray & table , TFile * afile = 0) ;

 void WriteHistogram(TH1* anhisto = 0, TString option = "");

 void WriteAllHistogramsInFile(TString filename = "result.root",
 			       TString option   = "recreate") ;
 
 void SaveHistogram(TH1* hist1, Char_t* title = "", Char_t* option = "", 
	            Bool_t logy = 0, Bool_t logx = 0, Int_t color = 5,
	            Int_t  stati = 111111, TString format = "ps");
		    
 void SetHistogramStyle(TString style = "seb");
 
 void SetSebStyle(float W = 0.3, float H = 0.15, 
 		  float titleW = 0.48, float titleH = 0.05);


 ClassDef(HistogramManager,1)// HistogramManager
};

#endif
