/**
 * MonitorOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 18, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include "Option/interface/MonitorOptionsDelegate.h"
#include "Option/interface/MiscOptions.h"
#include "Option/interface/MuonInJetOptions.h"
#include "Option/interface/PythiaOptions.h"
#include "Option/interface/TriggerOptions.h"

#include "Option/interface/MonitorOptions.h"

using s8::MonitorOptions;
using s8::MonitorOptionsDelegate;

MonitorOptions::MonitorOptions() throw()
{
    _delegate = 0;
}

MonitorOptions::~MonitorOptions() throw()
{
}

MonitorOptionsDelegate *MonitorOptions::delegate() const
{
    return _delegate;
}

void MonitorOptions::setDelegate(MonitorOptionsDelegate *delegate)
{
    _delegate = delegate;

    _misc_options->setDelegate(delegate);
    _muonInJetOptions->setDelegate(delegate);
    _pythiaOptions->setDelegate(delegate);
    _triggerOptions->setDelegate(delegate);
}

void MonitorOptions::init()
{
    _misc_options.reset(new MiscOptions());
    _misc_options->init();

    _muonInJetOptions.reset(new MuonInJetOptions());
    _muonInJetOptions->init();

    _pythiaOptions.reset(new PythiaOptions());
    _pythiaOptions->init();

    _triggerOptions.reset(new TriggerOptions());
    _triggerOptions->init();

    _description.reset(new po::options_description());
    _description->add(*_misc_options->description())
                 .add(*_muonInJetOptions->description())
                 .add(*_pythiaOptions->description())
                 .add(*_triggerOptions->description());
}

po::options_description *MonitorOptions::description() const
{
    return _description.get();
}

void MonitorOptions::print(std::ostream &out) const
{
    using std::endl;

    _misc_options->print(out);
    out << endl;

    _muonInJetOptions->print(out);
    out << endl;

    _pythiaOptions->print(out);
    out << endl;

    _triggerOptions->print(out);
}
