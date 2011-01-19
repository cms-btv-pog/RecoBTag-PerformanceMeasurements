/**
 * LeptonInJet
 * s8
 *
 * Created by Samvel Khalatian on Nov 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_LEPTON_IN_JET
#define S8_LEPTON_IN_JET

#include <vector>

namespace s8
{
    class Jet;
    class Lepton;
    class LeptonInJetDelegate;
    class PrimaryVertex;

    class LeptonInJet
    {
        public:
            typedef std::vector<Lepton *> Leptons;

            LeptonInJet() throw();
            virtual ~LeptonInJet() throw();

            LeptonInJetDelegate *delegate() const;
            void setDelegate(LeptonInJetDelegate *);

            void operator()(const PrimaryVertex *,
                            const Leptons *,
                            const Jet *);

        private:
            LeptonInJetDelegate *_delegate;
    };
}

#endif
