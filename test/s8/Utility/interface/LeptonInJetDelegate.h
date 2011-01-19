/**
 * LeptonInJetDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_LEPTON_IN_JET_DELEGATE
#define S8_LEPTON_IN_JET_DELEGATE

namespace s8
{
    class Jet;
    class Lepton;

    class LeptonInJetDelegate
    {
        public:
            virtual ~LeptonInJetDelegate() throw();

            virtual void leptonIsInJet(const Lepton *,
                                       const Jet *);

            virtual bool shouldLookForMoreLeptons();
    };
}

#endif
