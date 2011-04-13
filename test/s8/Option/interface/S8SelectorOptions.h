/**
 * S8SelectorOptions
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_S8SELECTOR_OPTIONS
#define S8_S8SELECTOR_OPTIONS

#include <string>

#include "FWCore/interface/Options.h"

namespace s8
{
    class S8SelectorOptionsDelegate;

    class S8SelectorOptions: public core::Options
    {
        public:
            S8SelectorOptions() throw();
            virtual ~S8SelectorOptions() throw();

            S8SelectorOptionsDelegate *delegate() const;
            void setDelegate(S8SelectorOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            void optionJetPtIsSet(const double &);
            void optionJetEtaIsSet(const std::string &);

            S8SelectorOptionsDelegate               *_delegate;
            std::auto_ptr<po::options_description>   _description;

            double _min_jet_pt;
            double _min_jet_eta;
            double _max_jet_eta;
    };
}

#endif
