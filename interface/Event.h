/**
 * Event 
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_EVENT
#define S8_EVENT

#include "RecoBTag/PerformanceMeasurements/interface/EventID.h"
#include "RecoBTag/PerformanceMeasurements/interface/Jet.h"
#include "RecoBTag/PerformanceMeasurements/interface/Muon.h"
#include "RecoBTag/PerformanceMeasurements/interface/PrimaryVertex.h"

namespace s8
{
    typedef std::vector<Jet>           JetCollection;
    typedef std::vector<Muon>          MuonCollection;
    typedef std::vector<PrimaryVertex> PrimaryVertexCollection;

    class Event
    {
        public:
            Event() throw();

            void reset();

            EventID &id();
            const EventID &id() const;

            JetCollection &jets();
            const JetCollection &jets() const;

            MuonCollection &muons();
            const MuonCollection &muons() const;

            PrimaryVertexCollection &primaryVertices();
            const PrimaryVertexCollection &primaryVertices() const;

        private:
            EventID _id;

            JetCollection           _jets;
            MuonCollection          _muons;
            PrimaryVertexCollection _primaryVertices;
    };

    inline Event::Event() throw()
    {
    }

    inline EventID &Event::id()
    {
        return _id;
    }

    inline const EventID &Event::id() const
    {
        return _id;
    }

    inline JetCollection &Event::jets()
    {
        return _jets;
    }

    inline const JetCollection &Event::jets() const
    {
        return _jets;
    }

    inline MuonCollection &Event::muons()
    {
        return _muons;
    }

    inline const MuonCollection &Event::muons() const
    {
        return _muons;
    }

    inline PrimaryVertexCollection &Event::primaryVertices()
    {
        return _primaryVertices;
    }

    inline const PrimaryVertexCollection &Event::primaryVertices() const
    {
        return _primaryVertices;
    }
}

#endif
