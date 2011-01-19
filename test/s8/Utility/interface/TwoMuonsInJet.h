/**
 * TwoMuonsInJet
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TWO_MUONS_IN_JET
#define S8_TWO_MUONS_IN_JET

#include <memory>

#include "Utility/interface/LeptonInJetDelegate.h"


namespace s8
{
    class Event;
    class Jet;
    class Lepton;
    class LeptonInJet;
    class TwoMuonsInJetDelegate;
    class TaggerOperatingPoint;


    class TwoMuonsInJet: public LeptonInJetDelegate
    {
        public:
            TwoMuonsInJet() throw();
            virtual ~TwoMuonsInJet() throw();

            TwoMuonsInJetDelegate *delegate() const;
            void setDelegate(TwoMuonsInJetDelegate *);

            void setAwayJetTaggerOperatingPoint(const TaggerOperatingPoint *);

            virtual void operator()(const Event *);

            // LeptonInJetDelegate interface
            //
            virtual void leptonIsInJet(const Lepton *,
                                       const Jet *);

            virtual bool shouldLookForMoreLeptons();

        private:
            TwoMuonsInJetDelegate *_delegate;

            std::auto_ptr<LeptonInJet>  _leptonInJet;
            const TaggerOperatingPoint *_awayJetTaggerOperatingPoint;

            const Lepton *_firstLepton;
            const Lepton *_secondLepton;
            const Jet *_jet;
    };
}

#endif
