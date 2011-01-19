/**
 * TwoMuonsInJetDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TWO_MUONS_IN_JET_DELEGATE
#define S8_TWO_MUONS_IN_JET_DELEGATE

namespace s8
{
    class Jet;
    class Lepton;

    // Always check the second muon pointer for ZERO value. There are events
    // when second muons can not be found in the jet
    //
    class TwoMuonsInJetDelegate
    {
        public:
            virtual ~TwoMuonsInJetDelegate() throw();

            virtual bool muonInJetShouldProcessJet(const Jet *);

            virtual void muonIsInJetPlusAwayJet(const Lepton *,
                                                const Lepton *,
                                                const Jet *);

            virtual void muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                      const Lepton *,
                                                      const Jet *);
    };
}

#endif
