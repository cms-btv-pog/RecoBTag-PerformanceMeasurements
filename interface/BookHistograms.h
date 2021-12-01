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
