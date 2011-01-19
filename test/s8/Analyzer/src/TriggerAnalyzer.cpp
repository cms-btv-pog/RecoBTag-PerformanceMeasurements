/**
 * TriggerAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <iostream>

#include "S8Tree/interface/S8TriggerCenter.h"

#include "Analyzer/interface/TriggerAnalyzer.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;

using s8::TriggerAnalyzer;

TriggerAnalyzer::TriggerAnalyzer() throw()
{
}

TriggerAnalyzer::~TriggerAnalyzer() throw()
{
}



// Analyzer interface
//
void TriggerAnalyzer::init()
{
}

void TriggerAnalyzer::treeDidLoad(const TreeInfo *,
                                  const TriggerCenter *triggerCenter)
{
    using std::setw;
    using std::left;

    cout << "Triggers List is saved in the log file" << endl
        << endl;

    clog << "Triggers:" << endl;
    const TriggerCenter::TriggerMap &triggers = triggerCenter->triggers();
    for(TriggerCenter::TriggerMap::const_iterator trigger = triggers.begin();
        triggers.end() != trigger;
        ++trigger)
    {
        clog << " > " << setw(20) << left << trigger->first << " "
            << trigger->second << endl;
    }
    clog << endl;
}

void TriggerAnalyzer::eventDidLoad(const Event *event)
{
}

void TriggerAnalyzer::save(TDirectory *directory) const
{
}
