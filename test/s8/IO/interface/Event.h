/**
 * Event
 * s8
 *
 * Created by Samvel Khalatian on Nov 21, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_EVENT
#define S8_EVENT

#include "S8Tree/interface/S8Fwd.h"

namespace s8
{
    class EventID;
    class GenEvent;

    class Event
    {
        public:
            Event() throw();

            const EventID *id() const;
            const GenEvent *gen() const;

            const Jets *jets() const;
            const Leptons *electrons() const;
            const Leptons *muons() const;
            const PrimaryVertices *primaryVertices() const;
            const Triggers *triggers() const;



            void setID(const EventID *);
            void setGen(const GenEvent *);
            void setJets(const Jets *);
            void setElectrons(const Leptons *);
            void setMuons(const Leptons *);
            void setPrimaryVertices(const PrimaryVertices *);
            void setTriggers(const Triggers *);

        private:
            const EventID *_id;
            const GenEvent *_gen;

            const Jets *_jets;
            const Leptons *_electrons;
            const Leptons *_muons;
            const PrimaryVertices *_primaryVertices;
            const Triggers *_triggers;
    };
};

#endif
