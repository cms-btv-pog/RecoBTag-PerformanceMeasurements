#include "TKey.h"
#include "TString.h"
#include "TFile.h"

void    MakeDirectory(TString pathnametofile );
TFile*  FindFile(TString pathnametofile );
Bool_t  FileExists(TString pathnametofile) ;
void    MoveFile(TFile* file ,TString oldname , TString pathname1 );
void    MoveFile(TString filename ,TString oldname , TString pathname1 );
void    CreateNewFile(TString path ,TString filename ,TString option = "" );


