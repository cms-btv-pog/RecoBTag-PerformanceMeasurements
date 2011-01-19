/**
 * PythiaOptions
 * s8
 *
 * Created by Samvel Khalatian on Dec 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PYTHIA_OPTIONS
#define S8_PYTHIA_OPTIONS

#include <string>

#include "FWCore/interface/Options.h"

namespace s8
{
    class PythiaOptionsDelegate;

    class PythiaOptions: public core::Options
    {
        public:
            PythiaOptions() throw();
            virtual ~PythiaOptions() throw();

            PythiaOptionsDelegate *delegate() const;
            void setDelegate(PythiaOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            void optionGluonSplittingIsSet(const std::string &);
            void optionPtHatIsSet(const std::string &);

            PythiaOptionsDelegate             *_delegate;
            std::auto_ptr<po::options_description>  _description;

            std::string _gluonSplitting;
            int         _minPtHat;
            int         _maxPtHat;
    };
}

#endif
