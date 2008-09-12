#ifndef Ptrel2DPlotter_h
#define Ptrel2DPlotter_h

#include "TH2F.h"

#include "BasePlotter.h"

class Ptrel2DPlotter : public BasePlotter {

 public:

  Ptrel2DPlotter(TString filename="") : BasePlotter(filename) {}

  void Write();

  ClassDef(Ptrel2DPlotter,1);

 protected:

  void Book(); 
  void Fill(BTagEvent*);

 private:

  TH1F * pTB_;
  TH1F * pTC_;
  TH1F * pTLight_;

  TH1F * deltaRB_;
  TH1F * deltaRC_;
  TH1F * deltaRLight_;

  TH1F * ntrkB_;
  TH1F * ntrkC_;
  TH1F * ntrkLight_;

  TH2F * pTrelVspTB_;
  TH2F * pTrelVspTC_;
  TH2F * pTrelVspTLight_;

  TH2F * pTrelVsDeltaRB_;
  TH2F * pTrelVsDeltaRC_;
  TH2F * pTrelVsDeltaRLight_;

  TH2F * pTrelVsNtrkB_;
  TH2F * pTrelVsNtrkC_;
  TH2F * pTrelVsNtrkLight_;

  TH2F * pTrelVsIps2trkB_;
  TH2F * pTrelVsIps2trkC_;
  TH2F * pTrelVsIps2trkLight_;

  TH2F * pTrelVsAbsIps2trkB_;
  TH2F * pTrelVsAbsIps2trkC_;
  TH2F * pTrelVsAbsIps2trkLight_;
 
};

#endif
