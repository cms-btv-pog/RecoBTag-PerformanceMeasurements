/**_________________________________________________________________
   class:   BTagEvent.cc
   package: RecoBTag/PerformanceMeasurements


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BTagEvent.cc,v 1.18 2008/09/12 16:41:30 bazterra Exp $

________________________________________________________________**/

#ifdef NOSCRAMV
#include "BTagEvent.h"
#else
#include "RecoBTag/PerformanceMeasurements/interface/BTagEvent.h"
#endif

#include <math.h>

ClassImp(BTagEvent)

// ROOT

//_______________________________________________________________
BTagEvent::BTagEvent()
{

    this->Reset();

}



//_______________________________________________________________
BTagEvent::~BTagEvent()
{
}

//_______________________________________________________________
void BTagEvent::Reset()
{

    event      = -1;
    run        = -1;
    evt_weight = 1.;
    // mc
    ptHat      = -1;
    // reco
    njets      = -1;
    nmuons     = -1;
    nelectrons = -1;

    nvertices  = -1;
    ngenjets   = -1;

    trackProbaVector_Size.clear();

    jet_pt.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_e.clear();
    jet_et.clear();
    jet_ntrks.clear();
    jet_flavour.clear();
    jet_hasLepton.clear();

    jetcorrection.clear();
    jet_Tracks_Probability.clear();


    genjet_pt.clear();
    genjet_eta.clear();
    genjet_phi.clear();
    genjet_e.clear();

    btag_TrkCounting_disc3D_1trk.clear();
    btag_TrkCounting_disc3D_2trk.clear();
    btag_TrkCounting_disc3D_3trk.clear();

    btag_TrkCounting_disc3D_1trk_is.clear();
    btag_TrkCounting_disc3D_2trk_is.clear();
    btag_TrkCounting_disc3D_3trk_is.clear();

    btag_ModTrkCounting_disc3D_1trk.clear();
    btag_ModTrkCounting_disc3D_2trk.clear();
    btag_ModTrkCounting_disc3D_3trk.clear();

    btag_ModTrkCounting_disc3D_1trk_is.clear();
    btag_ModTrkCounting_disc3D_2trk_is.clear();
    btag_ModTrkCounting_disc3D_3trk_is.clear();

    btag_JetProb_disc3D.clear();
    btag_negJetProb_disc3D.clear();
    btag_posJetProb_disc3D.clear();

    btag_NegTag_disc3D_1trk.clear();
    btag_NegTag_disc3D_2trk.clear();
    btag_NegTag_disc3D_3trk.clear();

    btag_NegTag_disc3D_1trk_is.clear();
    btag_NegTag_disc3D_2trk_is.clear();
    btag_NegTag_disc3D_3trk_is.clear();

    btag_SoftMuon_disc.clear();
    btag_SoftElectron_disc.clear();

    lepton.clear();
    tracks.clear();
}


std::vector< float > BTagEvent::getTrackProbabilies(std::vector< float > v, int ipType)
{

    std::vector< float > vectTrackProba;

    for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); q++)
    {
        //positives and negatives tracks
        double p3d = 0;
        if (ipType == 0)
        {
            if ( *q>=0)
            {
                p3d= (*q)/2.;
            }
            else
            {
                p3d=1.+(*q)/2.;
            }
            //if(-log(p3d)> 5) p3d=exp(-5.0);
            vectTrackProba.push_back(p3d);
        }

        //positives tracks only
        if (ipType == 1 && *q >=0 )
        {
            vectTrackProba.push_back(*q);
        }

        //negatives tracks only
        if (ipType == 2 && *q <0)
        {
            vectTrackProba.push_back(-(*q));
        }




    }


    return vectTrackProba;


}






double BTagEvent::calculProbability(std::vector< float > v)
{

    int ngoodtracks=v.size();
    double SumJet=0.;
    double m_minTrackProb = 0.005;
    for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); q++)
    {
        SumJet+=(*q>m_minTrackProb)?log(*q):log(m_minTrackProb);
    }

    double ProbJet;
    double Loginvlog=0;

    if (SumJet<0.)
    {
        if (ngoodtracks>=2)
        {
            Loginvlog=log(-SumJet);
        }
        double Prob=1.;
        double lfact=1.;
        for (int l=1; l!=ngoodtracks; l++)
        {
            lfact*=l;
            Prob+=exp(l*Loginvlog-log(1.*lfact));
        }
        double LogProb=log(Prob);
        ProbJet=
            std::min(exp(std::max(LogProb+SumJet,-30.)),1.);
    }
    else
    {
        ProbJet=1.;
    }
    //if(ProbJet>1)
    //std::cout << "ProbJet too high: "  << ProbJet << std::endl;

    //double LogProbJet=-log(ProbJet);
    //  //return 1.-ProbJet;
    return -log10(ProbJet)/4.;

}

