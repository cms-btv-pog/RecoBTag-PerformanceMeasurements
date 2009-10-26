#ifndef BTagEvent_h
#define BTagEvent_h

/**_________________________________________________________________
   class:   BTagEvent.h
   package: RecoBTag/PerformanceMeasurements


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagEvent.h,v 1.18 2008/09/12 16:41:30 bazterra Exp $

________________________________________________________________**/


// C++ standard

#include <vector>

// ROOT
#include "TObject.h"
#include "TMatrixDSym.h"

#ifdef NOSCRAMV
#include "BTagTrackEvent.h"
#include "BTagLeptonEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagTrackEvent.h"
#include "RecoBTag/PerformanceMeasurements/interface/BTagLeptonEvent.h"
#endif


class BTagEvent : public TObject
{

public:

    BTagEvent();
    ~BTagEvent();

    void                 Reset();
    double               calculProbability(std::vector< float > );
    std::vector< float > getTrackProbabilies(std::vector< float > v , int ipType);

    Int_t event;       // event number
    Int_t run;         // run number
    double evt_weight; // event weight
    float ptHat;       // hatted pt if runnging over the soup

    // reco
    Int_t njets;       // number of jets
    Int_t nmuons;      // number of muons
    Int_t nelectrons;  // number of electrons
    Int_t nvertices;   // number of vertices
    Int_t ngenjets;    // number of generated jets
    std::vector< Int_t > trackProbaVector_Size;  // size of the vector of tracks
    //== track multiplicity of quality cuted tracks


    std::vector< float > jet_pt;
    std::vector< float > jet_eta;
    std::vector< float > jet_phi;
    std::vector< float > jet_e;
    std::vector< float > jet_et;
    std::vector< int > jet_ntrks;

    //std::vector< int > jet_flavour_phy;
    std::vector< int > jet_flavour;
    std::vector< int > jet_hasLepton;
    std::vector< float > jetcorrection;
    std::vector< std::vector< float > > jet_Tracks_Probability;
    std::vector< float > genjet_pt;
    std::vector< float > genjet_eta;
    std::vector< float > genjet_phi;
    std::vector< float > genjet_e;

    std::vector< float > btag_TrkCounting_disc3D_1trk;
    std::vector< float > btag_TrkCounting_disc3D_2trk;
    std::vector< float > btag_TrkCounting_disc3D_3trk;

    std::vector<std::vector<bool> > btag_TrkCounting_disc3D_1trk_is;
    std::vector<std::vector<bool> > btag_TrkCounting_disc3D_2trk_is;
    std::vector<std::vector<bool> > btag_TrkCounting_disc3D_3trk_is;

    std::vector< float > btag_ModTrkCounting_disc3D_1trk;
    std::vector< float > btag_ModTrkCounting_disc3D_2trk;
    std::vector< float > btag_ModTrkCounting_disc3D_3trk;

    std::vector<std::vector<bool> > btag_ModTrkCounting_disc3D_1trk_is;
    std::vector<std::vector<bool> > btag_ModTrkCounting_disc3D_2trk_is;
    std::vector<std::vector<bool> > btag_ModTrkCounting_disc3D_3trk_is;

    std::vector< float > btag_JetProb_disc3D;
    std::vector< float > btag_negJetProb_disc3D;
    std::vector< float > btag_posJetProb_disc3D;

    std::vector< float > btag_NegTag_disc3D_1trk;
    std::vector< float > btag_NegTag_disc3D_2trk;
    std::vector< float > btag_NegTag_disc3D_3trk;

    std::vector<std::vector<bool> > btag_NegTag_disc3D_1trk_is;
    std::vector<std::vector<bool> > btag_NegTag_disc3D_2trk_is;
    std::vector<std::vector<bool> > btag_NegTag_disc3D_3trk_is;

    std::vector< float > btag_SoftMuon_disc;
    std::vector< float > btag_SoftElectron_disc;

    std::vector< BTagTrackEvent > tracks;
    std::vector< BTagLeptonEvent > lepton;

    ClassDef(BTagEvent,1);

};

#endif
