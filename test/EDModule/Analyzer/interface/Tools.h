/**
 * Tools
 * top
 *
 * Created by Samvel Khalatian on Sep 7, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef TOP_EDANALYZER_TOOLS
#define TOP_EDANALYZER_TOOLS

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "Tree/Top/interface/TopElectron.h"
#include "Tree/Top/interface/TopMuon.h"

class TLorentzVector;
class TVector3;

namespace reco
{
    class MuonIsolation;
}

namespace top
{
    class Jet;
    class JetEnergy;
    class ElectronIsolation;
    class MuonIsolation;

    namespace tools
    {
        void setP4(TLorentzVector *,
                   const math::XYZTLorentzVector *);

        inline void setP4(TLorentzVector *p4_1,
                          const math::XYZTLorentzVector &p4_2)
        {
            setP4(p4_1, &p4_2);
        }

        inline void setP4(TLorentzVector &p4_1,
                          const math::XYZTLorentzVector &p4_2)
        {
            setP4(&p4_1, &p4_2);
        }

        void setVertex(TVector3 &,
                       const math::XYZPoint &);

        inline void setVertex(TVector3 *v1,
                              const math::XYZPoint &v2)
        {
            setVertex(*v1, v2);
        }

        void setEnergy(top::Jet &, const reco::CaloJet::Specific &);

        void setIsolation(top::Muon &,
                          const top::Muon::ISO &,
                          const reco::MuonIsolation &);

        void setIsolation(top::Electron &,
                          const top::Electron::ISO &,
                          const reco::GsfElectron::IsolationVariables &);
    }
}

#endif
