/**
 * PtRelSelectorOptions
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PTREL_SELECTOR_OPTIONS
#define S8_PTREL_SELECTOR_OPTIONS

#include <string>

#include "FWCore/interface/Options.h"

namespace s8
{
    class PtRelSelectorOptionsDelegate;

    class PtRelSelectorOptions: public core::Options
    {
        public:
            PtRelSelectorOptions() throw();
            virtual ~PtRelSelectorOptions() throw();

            PtRelSelectorOptionsDelegate *delegate() const;
            void setDelegate(PtRelSelectorOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            void optionMuonPtIsSet(const double &);

            PtRelSelectorOptionsDelegate           *_delegate;
            std::auto_ptr<po::options_description>  _description;

            double      _muonPt;
    };
}

#endif
