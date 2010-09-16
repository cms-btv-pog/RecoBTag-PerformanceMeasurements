//
//
// Package:    RecoBTag/PerformanceMeasurements
// Class:      PM_Histograms
//
/**\class PerformanceMeasurements/PM_Histograms

 Description:

         Author: Francisco Yumiceva, Fermilab
*/
//
// $Id: PMHistograms.h,v 1.5 2010/09/14 20:36:04 samvel Exp $
//
//


#ifndef PMHistograms_H
#define PMHistograms_H

#include "RecoBTag/PerformanceMeasurements/interface/TH1Store.h"
#include "TH2F.h"
#include "TString.h"
#include "TLorentzVector.h"

class PMHistograms
{
public:

    //PMHistograms() {}
    PMHistograms( TH1Store *hstore )
    {
        fstore = hstore;
    }
    ~PMHistograms() { }

    void Add();
    void FillHistos(const std::string &type,
                    const TLorentzVector &p4MuJet,
                    const double &ptrel,
                    int JetFlavor, bool tagged=false);

private:
    void FillHisto(const std::string &flavor,
                   const std::string &type,
                   const TLorentzVector &p4,
                   const double &ptrel,
                   const bool &tagged);

    void Fill(const std::string &prefix,
              const std::string &suffix,
              const TLorentzVector &p4,
              const double &ptrel);

    TH1Store *fstore;

};

#endif

