/**
 * MonitorOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 19, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MONITOR_OPTIONS_DELEGATE
#define S8_MONITOR_OPTIONS_DELEGATE

#include "Option/interface/MuonInJetOptionsDelegate.h"
#include "Option/interface/PythiaOptionsDelegate.h"
#include "Option/interface/TriggerOptionsDelegate.h"

namespace s8
{
    class MonitorOptionsDelegate: public MuonInJetOptionsDelegate,
                                  public PythiaOptionsDelegate,
                                  public TriggerOptionsDelegate
    {
        public:
            MonitorOptionsDelegate() throw();
            virtual ~MonitorOptionsDelegate() throw();
    };
}

#endif
