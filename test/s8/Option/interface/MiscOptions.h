/**
 * MiscOptions
 * s8
 *
 * Created by Samvel Khalatian on Mar 9, 2011
 * Copyright 2011, All rights reserved
 */

#ifndef S8_MICS_OPTIONS
#define S8_MICS_OPTIONS

#include <string>

#include "FWCore/interface/Options.h"
#include "Utility/interface/Range.h"

namespace s8
{
    class MiscOptionsDelegate;

    class MiscOptions: public core::Options
    {
        public:
            MiscOptions() throw();
            virtual ~MiscOptions() throw();

            MiscOptionsDelegate *delegate() const;
            void setDelegate(MiscOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            void optionPrimaryVerticesIsSet(const std::string &);

            MiscOptionsDelegate *_delegate;
            std::auto_ptr<po::options_description> _description;

            Range _primary_vertices;
    };
}

#endif
