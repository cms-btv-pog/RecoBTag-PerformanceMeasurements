#ifndef PFTools_h
#define PFTools_h

/**_________________________________________________________________
   class:   PFTools.h

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: PFTools.h,v 1.6 2010/03/09 16:50:43 jindal Exp $

________________________________________________________________**/

#include <iostream>
#include <map>
#include <string>

//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"


class WorkingPoint
{
public:

    WorkingPoint() {}

    WorkingPoint(
        edm::InputTag t, std::string n, double min, double max = 0
    ) : intag_(t), alias_(n), min_(min), max_(max) {}

    WorkingPoint(
        edm::InputTag t, std::string n, double min, double max , std::map<std::string, double> const & list
    ) : intag_(t), alias_(n), min_(min), max_(max)
    {
        if ( !list.empty() ) wpmap_ = list;
    }

    std::map<std::string, double > const & list() const
    {
        return wpmap_;
    }

    edm::InputTag inputTag() const
    {
        return intag_;
    }

    std::string alias () const
    {
        return alias_;
    }

    //std::string name () const
    //{
    //    return alias_;
    //}

    double cut() const
    {
        return min_;
    }

    double Minimum() const
    {
        return min_;
    }

    double Maximum() const
    {
        return max_;
    }

    void print () const
    {
        std::cout << " Working point: collection: " << inputTag() << " alias: " << alias() << " discriminator min: " << Minimum() << " max: " << Maximum() << std::endl;
        for (std::map<std::string, double >::const_iterator icut = wpmap_.begin(); icut != wpmap_.end(); ++icut)
        {
            std::string aalias = icut->first;
            double cut = icut->second;
            std::cout << "       name: " << aalias << " cut: " << cut << std::endl;
        }
    }

private:

    edm::InputTag intag_;
    std::string alias_;
    double min_;
    double max_;
    std::map< std::string, double > wpmap_;

};


#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "FWCore/Framework/interface/Event.h"

namespace PFTools
{

reco::GenJet GetGenJet(reco::CaloJet, reco::GenJetCollection);

int TaggedJet(reco::CaloJet const &, edm::Handle<reco::JetTagCollection> const &);

int TaggedJet(reco::CaloJet const &, edm::Handle<std::vector<reco::TrackIPTagInfo> > const &);

std::map<std::string, bool> GetBTaggingMap(std::vector<WorkingPoint> const &, reco::CaloJet const &, edm::Event const &, double ptrel = 0);

// THis function seems to not be use, if that the case please remove it !!!
// SimTrack GetGenTrk(reco::Track const &, edm::SimTrackContainer const *, edm::SimVertexContainer const *);

}

#endif
