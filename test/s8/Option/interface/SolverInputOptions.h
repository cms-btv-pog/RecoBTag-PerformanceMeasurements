/**
 * SolverInputOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_SOLVER_INPUT_OPTIONS
#define S8_SOLVER_INPUT_OPTIONS

#include <string>

#include "FWCore/interface/Options.h"

namespace s8
{
    class MuonInJetOptions;
    class MiscOptions;
    class PythiaOptions;
    class SolverInputOptionsDelegate;
    class TriggerOptions;

    class SolverInputOptions: public core::Options
    {
        public:
            SolverInputOptions() throw();
            virtual ~SolverInputOptions() throw();

            SolverInputOptionsDelegate *delegate() const;
            void setDelegate(SolverInputOptionsDelegate *);

            // Options interface
            //
            virtual void init();

            virtual po::options_description *description() const;

            virtual void print(std::ostream &out) const;

        private:
            void optionDataIsSet(const bool &);

            SolverInputOptionsDelegate             *_delegate;

            std::auto_ptr<po::options_description>  _description;
            std::auto_ptr<po::options_description>  _hiddenDescription;

            std::auto_ptr<MiscOptions> _misc_options;
            std::auto_ptr<MuonInJetOptions> _muonInJetOptions;
            std::auto_ptr<PythiaOptions>    _pythiaOptions;
            std::auto_ptr<TriggerOptions>   _triggerOptions;

            bool        _isData;
    };
}

#endif
