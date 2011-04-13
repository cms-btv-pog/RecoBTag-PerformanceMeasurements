/**
 * MuonInJetOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MUON_IN_JET_OPTIONS
#define S8_MUON_IN_JET_OPTIONS

#include <string>

#include "FWCore/interface/Options.h"
#include "Utility/interface/Range.h"

namespace s8
{
    class MuonInJetOptionsDelegate;

    class MuonInJetOptions: public core::Options
    {
        public:
            MuonInJetOptions() throw();
            virtual ~MuonInJetOptions() throw();

            MuonInJetOptionsDelegate *delegate() const;
            void setDelegate(MuonInJetOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            void optionTagIsSet(const std::string &);
            void optionAwayTagIsSet(const std::string &);
            void optionMuonPtIsSet(const std::string &);
            void optionJetPtIsSet(const std::string &);
            void optionJetEtaIsSet(const std::string &);

            MuonInJetOptionsDelegate               *_delegate;
            std::auto_ptr<po::options_description>  _description;

            std::string _tag;
            std::string _awayTag;
            Range       _muon_pt;
            Range       _jet_pt;
            Range       _jet_eta;
    };
}

#endif
