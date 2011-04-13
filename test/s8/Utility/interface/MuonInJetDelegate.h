/**
 * MuonInJetDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MUON_IN_JET_DELEGATE
#define S8_MUON_IN_JET_DELEGATE

namespace s8
{
    class Jet;
    class Lepton;

    class MuonInJetDelegate
    {
        public:
            virtual ~MuonInJetDelegate() throw();

            virtual bool muonInJetShouldProcessJet(const Jet *);

            virtual bool shouldSkipMuonInJetPlusAwayJet(const Lepton *,
                                                        const Jet *);

            virtual void muonIsInJetPlusAwayJet(const Lepton *,
                                                const Jet *);

            virtual void muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                      const Jet *);
    };
}

#endif
