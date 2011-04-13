/**
 * MuonInJet
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MUON_IN_JET
#define S8_MUON_IN_JET

#include <memory>

#include "Utility/interface/LeptonInJetDelegate.h"

namespace s8
{
    class Event;
    class Jet;
    class Lepton;
    class LeptonInJet;
    class MuonInJetDelegate;
    class TaggerOperatingPoint;

    class MuonInJet: public LeptonInJetDelegate
    {
        public:
            MuonInJet() throw();
            virtual ~MuonInJet() throw();

            MuonInJetDelegate *delegate() const;
            void setDelegate(MuonInJetDelegate *);

            void setAwayJetTaggerOperatingPoint(const TaggerOperatingPoint *);

            virtual void operator()(const Event *);

            // LeptonInJetDelegate interface
            //
            virtual void leptonIsInJet(const Lepton *,
                                       const Jet *);

            virtual bool shouldLookForMoreLeptons();

        private:
            MuonInJetDelegate *_delegate;

            std::auto_ptr<LeptonInJet>  _leptonInJet;
            const TaggerOperatingPoint *_awayJetTaggerOperatingPoint;

            const Lepton *_lepton;
            const Jet    *_jet;
    };
}

#endif
