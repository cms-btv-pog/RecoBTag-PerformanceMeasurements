/**
 * Muon
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MUON
#define S8_MUON

#include <utility>

#include <TLorentzVector.h>
#include <TVector3.h>

#include "RecoBTag/PerformanceMeasurements/interface/GenParticle.h"

namespace s8
{
    class Muon
    {
        public:
            // Value and associated Error pair
            //
            typedef std::pair<double,double> ImpactParameter;

            Muon() throw();

            ImpactParameter &impactParameter();
            const ImpactParameter &impactParameter() const;

            TLorentzVector &p4();
            const TLorentzVector &p4() const;

            TVector3 &vertex();
            const TVector3 &vertex() const;

            GenParticle &genParticle();
            const GenParticle &genParticle() const;

        private:
            ImpactParameter _impactParameter;
            TLorentzVector  _p4;
            TVector3        _vertex;

            GenParticle     _genParticle;
    };

    inline Muon::Muon() throw()
    {
    }

    inline Muon::ImpactParameter &Muon::impactParameter()
    {
        return _impactParameter;
    }

    inline const Muon::ImpactParameter &Muon::impactParameter() const
    {
        return _impactParameter;
    }

    inline TLorentzVector &Muon::p4()
    {
        return _p4;
    }

    inline const TLorentzVector &Muon::p4() const
    {
        return _p4;
    }

    inline TVector3 &Muon::vertex()
    {
        return _vertex;
    }

    inline const TVector3 &Muon::vertex() const
    {
        return _vertex;
    }

    inline GenParticle &Muon::genParticle()
    {
        return _genParticle;
    }

    inline const GenParticle &Muon::genParticle() const
    {
        return _genParticle;
    }
}

#endif
