/**
 * SolverInputOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_SOLVER_INPUT_OPTIONS_DELEGATE
#define S8_SOLVER_INPUT_OPTIONS_DELEGATE

#include "Option/interface/MiscOptionsDelegate.h"
#include "Option/interface/MuonInJetOptionsDelegate.h"
#include "Option/interface/PythiaOptionsDelegate.h"
#include "Option/interface/TriggerOptionsDelegate.h"

namespace s8
{
    class SolverInputOptionsDelegate: public MiscOptionsDelegate,
                                      public MuonInJetOptionsDelegate,
                                      public PythiaOptionsDelegate,
                                      public TriggerOptionsDelegate
    {
        public:
            SolverInputOptionsDelegate() throw();
            virtual ~SolverInputOptionsDelegate() throw();

            virtual void optionDataIsSet(const bool &);
    };
}

#endif
