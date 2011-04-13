/**
 * MuonInJetOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MUON_IN_JET_OPTIONS_DELEGATE
#define S8_MUON_IN_JET_OPTIONS_DELEGATE

#include <string>

namespace s8
{
    class Range;

    class MuonInJetOptionsDelegate
    {
        public:
            MuonInJetOptionsDelegate() throw();
            virtual ~MuonInJetOptionsDelegate() throw();

            virtual void optionTagIsSet(const std::string &);
            virtual void optionAwayTagIsSet(const std::string &);
            virtual void optionMuonPtIsSet(const Range &);
            virtual void optionJetPtIsSet(const Range &);
            virtual void optionJetEtaIsSet(const Range &);
    };
}

#endif
