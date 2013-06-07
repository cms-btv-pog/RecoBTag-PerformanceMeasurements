#ifndef BOOKHISTOGRAMS_H
#define BOOKHISTOGRAMS_H

#include <TH1F.h>

#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class BookHistograms {

  public :

    BookHistograms (TFileDirectory dir) ;
    ~BookHistograms () {} ;

    TH1F* hData_All_NJets       ;
    TH1F* hData_All_NTracks     ;
    TH1F* hData_All_JetPt       ;
    TH1F* hData_All_JetEta      ;
    TH1F* hData_NJets           ;
    TH1F* hData_NTracks         ;
    TH1F* hData_JetPt           ;
    TH1F* hData_JetEta          ;
    TH1F* hData_Tagger          ;
    TH1F* hData_Tagger_TCHE     ;
    TH1F* hData_Tagger_TCHP     ;
    TH1F* hData_Tagger_JP	      ;
    TH1F* hData_Tagger_SSVHE    ;
    TH1F* hData_Tagger_SSVHP    ;
    TH1F* hData_Tagger_CSV      ;
    TH1F* hData_Tagger_MU	;

    TH1F* hAllFlav_Flavour      ;
    TH1F* hAllFlav_Tagger       ;
    TH1F* hAllFlav_Tagger_Gam   ;
    TH1F* hAllFlav_Tagger_K0s   ;
    TH1F* hAllFlav_Tagger_Lam   ;
    TH1F* hAllFlav_Tagger_Bwd   ;
    TH1F* hAllFlav_Tagger_Cwd   ;
    TH1F* hAllFlav_Tagger_Tau   ;
    TH1F* hAllFlav_Tagger_Int   ;
    TH1F* hAllFlav_Tagger_Fak   ;
    TH1F* hAllFlav_Tagger_Bad   ;
    TH1F* hAllFlav_Tagger_Oth   ;

    TH1F* hLightFlav_Tagger     ;
    TH1F* hGluonFlav_Tagger     ;
    TH1F* hUDSFlav_Tagger       ;
    TH1F* hCFlav_Tagger         ;
    TH1F* hBFlav_Tagger         ;

    TH1F*  IPSign_cat0          ;
    TH1F*  IPSign_cat1          ;
    TH1F*  IPSign_cat2          ;
    TH1F*  IPSign_cat3          ;
    TH1F*  IPSign_cat4          ;
    TH1F*  IPSign_cat5          ;
    TH1F*  IPSign_cat6          ;
    TH1F*  IPSign_cat7          ;
    TH1F*  IPSign_cat8          ;
    TH1F*  IPSign_cat9          ;
    TH1F*  TrackProbaNeg        ;
    TH1F*  TrackProbaNeg_Cat0   ;
    TH1F*  TrackProbaNeg_Cat1   ;
    TH1F*  TrackProbaNeg_Cat2   ;
    TH1F*  TrackProbaNeg_Cat3   ;
    TH1F*  TrackProbaNeg_Cat4   ;
    TH1F*  TrackProbaNeg_Cat5   ;
    TH1F*  TrackProbaNeg_Cat6   ;
    TH1F*  TrackProbaNeg_Cat7   ;
    TH1F*  TrackProbaNeg_Cat8   ;
    TH1F*  TrackProbaNeg_Cat9   ;
    TH1F*  TrackProbJet80       ;
    TH1F*  TrackProbJet80_Cat0  ;
    TH1F*  TrackProbJet80_Cat1  ;
    TH1F*  TrackProbJet80_Cat2  ;
    TH1F*  TrackProbJet80_Cat3  ;
    TH1F*  TrackProbJet80_Cat4  ;
    TH1F*  TrackProbJet80_Cat5  ;
    TH1F*  TrackProbJet80_Cat6  ;
    TH1F*  TrackProbJet80_Cat7  ;
    TH1F*  TrackProbJet80_Cat8  ;
    TH1F*  TrackProbJet80_Cat9  ;

};

#endif

#ifdef BOOKHISTOGRAMS_H

BookHistograms::BookHistograms (TFileDirectory dir) {

  ///////////////
  // Some Histograms

  hData_All_NJets       = dir.make<TH1F>("hData_All_NJets","nb. of jets"     ,21, -0.5, 20.5);
  hData_All_NTracks     = dir.make<TH1F>("hData_All_NTracks","nb. of tracks" ,20,  0.5, 25.5);
  hData_All_JetPt       = dir.make<TH1F>("hData_All_JetPt","pt(jet)"         ,50,  20., 520.);
  hData_All_JetEta      = dir.make<TH1F>("hData_All_JetEta","|#eta(jet)|"    ,24,   0., 2.5 );
  hData_NJets           = dir.make<TH1F>("hData_NJets","nb. of jets"         ,21, -0.5, 20.5);
  hData_NTracks         = dir.make<TH1F>("hData_NTracks","nb. of tracks"     ,20,  0.5, 25.5);
  hData_JetPt           = dir.make<TH1F>("hData_JetPt","pt(jet)"             ,50,  20., 520.);
  hData_JetEta          = dir.make<TH1F>("hData_JetEta","|#eta(jet)|"        ,24,  0. , 2.5 );

  hData_Tagger  	      = dir.make<TH1F>("hData_Tagger","Tagger"             ,100,-25.,25.);
  hData_Tagger_TCHE	    = dir.make<TH1F>("hData_Tagger_TCHE","Tagger_TCHE"   ,100,-25.,25.);
  hData_Tagger_TCHP	    = dir.make<TH1F>("hData_Tagger_TCHP","Tagger_TCHP"   ,100,-25.,25.);
  hData_Tagger_JP	      = dir.make<TH1F>("hData_Tagger_JP","Tagger_JP"       ,100,-25.,25.);
  hData_Tagger_SSVHE    = dir.make<TH1F>("hData_Tagger_SSVHE ","Tagger_SSVHE",100,-25.,25.);
  hData_Tagger_SSVHP    = dir.make<TH1F>("hData_Tagger_SSVHP","Tagger_SSVHP" ,100,-25.,25.);
  hData_Tagger_CSV	    = dir.make<TH1F>("hData_Tagger_CSV","Tagger_CSV"     ,100,-25.,25.);
  hData_Tagger_MU	      = dir.make<TH1F>("hData_Tagger_MU","Tagger_MU"       ,100,-25.,25.);

  hAllFlav_Flavour      = dir.make<TH1F>("hAllFlav_Flavour","Flavour"   ,22,-0.5,21.5);
  hAllFlav_Tagger	      = dir.make<TH1F>("hAllFlav_Tagger","Tagger"     ,100,-25.,25.);
  hAllFlav_Tagger_Gam	  = dir.make<TH1F>("hAllFlav_Tagger_Gam","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_K0s	  = dir.make<TH1F>("hAllFlav_Tagger_K0s","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Lam	  = dir.make<TH1F>("hAllFlav_Tagger_Lam","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Bwd	  = dir.make<TH1F>("hAllFlav_Tagger_Bwd","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Cwd	  = dir.make<TH1F>("hAllFlav_Tagger_Cwd","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Tau	  = dir.make<TH1F>("hAllFlav_Tagger_Tau","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Int	  = dir.make<TH1F>("hAllFlav_Tagger_Int","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Fak	  = dir.make<TH1F>("hAllFlav_Tagger_Fak","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Bad	  = dir.make<TH1F>("hAllFlav_Tagger_Bad","Tagger" ,100,-25.,25.);
  hAllFlav_Tagger_Oth	  = dir.make<TH1F>("hAllFlav_Tagger_Oth","Tagger" ,100,-25.,25.);

  hLightFlav_Tagger	    = dir.make<TH1F>("hLightFlav_Tagger","Tagger",100,-25.,25.);
  hGluonFlav_Tagger	    = dir.make<TH1F>("hGluonFlav_Tagger","Tagger",100,-25.,25.);
  hUDSFlav_Tagger	      = dir.make<TH1F>("hUDSFlav_Tagger","Tagger"  ,100,-25.,25.);
  hCFlav_Tagger 	      = dir.make<TH1F>("hCFlav_Tagger","Tagger"    ,100,-25.,25.);
  hBFlav_Tagger 	      = dir.make<TH1F>("hBFlav_Tagger","Tagger"    ,100,-25.,25.);

  IPSign_cat0		        = dir.make<TH1F>("IPSign_cat0","-IP/#sigma",50, 0., 25.);
  IPSign_cat1		        = dir.make<TH1F>("IPSign_cat1","-IP/#sigma",50, 0., 25.);
  IPSign_cat2		        = dir.make<TH1F>("IPSign_cat2","-IP/#sigma",50, 0., 25.);
  IPSign_cat3		        = dir.make<TH1F>("IPSign_cat3","-IP/#sigma",50, 0., 25.);
  IPSign_cat4		        = dir.make<TH1F>("IPSign_cat4","-IP/#sigma",50, 0., 25.);
  IPSign_cat5		        = dir.make<TH1F>("IPSign_cat5","-IP/#sigma",50, 0., 25.);
  IPSign_cat6		        = dir.make<TH1F>("IPSign_cat6","-IP/#sigma",50, 0., 25.);
  IPSign_cat7		        = dir.make<TH1F>("IPSign_cat7","-IP/#sigma",50, 0., 25.);
  IPSign_cat8		        = dir.make<TH1F>("IPSign_cat8","-IP/#sigma",50, 0., 25.);
  IPSign_cat9		        = dir.make<TH1F>("IPSign_cat9","-IP/#sigma",50, 0., 25.);

  TrackProbaNeg       	= dir.make<TH1F>("TrackProbaNEG","TrackProbaNEG",51, 0., 1.02);
  TrackProbaNeg_Cat0	  = dir.make<TH1F>("TrackProbaNEG Cat0","TrackProbaNEG Cat0",51, 0., 1.02);
  TrackProbaNeg_Cat1	  = dir.make<TH1F>("TrackProbaNEG Cat1","TrackProbaNEG Cat1",51, 0., 1.02);
  TrackProbaNeg_Cat2	  = dir.make<TH1F>("TrackProbaNEG Cat2","TrackProbaNEG Cat2",51, 0., 1.02);
  TrackProbaNeg_Cat3	  = dir.make<TH1F>("TrackProbaNEG Cat3","TrackProbaNEG Cat3",51, 0., 1.02);
  TrackProbaNeg_Cat4	  = dir.make<TH1F>("TrackProbaNEG Cat4","TrackProbaNEG Cat4",51, 0., 1.02);
  TrackProbaNeg_Cat5	  = dir.make<TH1F>("TrackProbaNEG Cat5","TrackProbaNEG Cat5",51, 0., 1.02);
  TrackProbaNeg_Cat6	  = dir.make<TH1F>("TrackProbaNEG Cat6","TrackProbaNEG Cat6",51, 0., 1.02);
  TrackProbaNeg_Cat7	  = dir.make<TH1F>("TrackProbaNEG Cat7","TrackProbaNEG Cat7",51, 0., 1.02);
  TrackProbaNeg_Cat8	  = dir.make<TH1F>("TrackProbaNEG Cat8","TrackProbaNEG Cat8",51, 0., 1.02);
  TrackProbaNeg_Cat9	  = dir.make<TH1F>("TrackProbaNEG Cat9","TrackProbaNEG Cat9",51, 0., 1.02);
  TrackProbJet80	      = dir.make<TH1F>("TrackProbJet80","TrackProbJet80",51, 0., 1.02);
  TrackProbJet80_Cat0	  = dir.make<TH1F>("TrackProbJet80 Cat0","TrackProbJet80 Cat0",51, 0., 1.02);
  TrackProbJet80_Cat1	  = dir.make<TH1F>("TrackProbJet80 Cat1","TrackProbJet80 Cat1",51, 0., 1.02);
  TrackProbJet80_Cat2	  = dir.make<TH1F>("TrackProbJet80 Cat2","TrackProbJet80 Cat2",51, 0., 1.02);
  TrackProbJet80_Cat3	  = dir.make<TH1F>("TrackProbJet80 Cat3","TrackProbJet80 Cat3",51, 0., 1.02);
  TrackProbJet80_Cat4	  = dir.make<TH1F>("TrackProbJet80 Cat4","TrackProbJet80 Cat4",51, 0., 1.02);
  TrackProbJet80_Cat5	  = dir.make<TH1F>("TrackProbJet80 Cat5","TrackProbJet80 Cat5",51, 0., 1.02);
  TrackProbJet80_Cat6	  = dir.make<TH1F>("TrackProbJet80 Cat6","TrackProbJet80 Cat6",51, 0., 1.02);
  TrackProbJet80_Cat7	  = dir.make<TH1F>("TrackProbJet80 Cat7","TrackProbJet80 Cat7",51, 0., 1.02);
  TrackProbJet80_Cat8	  = dir.make<TH1F>("TrackProbJet80 Cat8","TrackProbJet80 Cat8",51, 0., 1.02);
  TrackProbJet80_Cat9	  = dir.make<TH1F>("TrackProbJet80 Cat9","TrackProbJet80 Cat9",51, 0., 1.02);

}

#endif
