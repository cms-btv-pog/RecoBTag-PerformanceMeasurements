/**
 * SolverInputOptions
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#include <iomanip>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "Option/interface/MuonInJetOptions.h"
#include "Option/interface/PythiaOptions.h"
#include "Option/interface/SolverInputOptionsDelegate.h"

#include "Option/interface/SolverInputOptions.h"

using std::runtime_error;
using std::string;

using boost::lexical_cast;
using boost::regex;
using boost::smatch;

using s8::SolverInputOptions;
using s8::SolverInputOptionsDelegate;

SolverInputOptions::SolverInputOptions() throw():
    _isData(false)
{
    _delegate = 0;
}

SolverInputOptions::~SolverInputOptions() throw()
{
}

SolverInputOptionsDelegate *SolverInputOptions::delegate() const
{
    return _delegate;
}

void SolverInputOptions::setDelegate(SolverInputOptionsDelegate *delegate)
{
    _delegate = delegate;

    _muonInJetOptions->setDelegate(delegate);
    _pythiaOptions->setDelegate(delegate);
}

void SolverInputOptions::init()
{
    _hiddenDescription.reset(new po::options_description("Solver Input Options"));
    _hiddenDescription->add_options()
        ("data",
            po::value<bool>()->default_value(false)->notifier(
                boost::bind(&SolverInputOptions::optionDataIsSet,
                            this, _1)),
            "Input is Data if set (no flavour splitting)")
    ;

    _muonInJetOptions.reset(new MuonInJetOptions());
    _muonInJetOptions->init();

    _pythiaOptions.reset(new PythiaOptions());
    _pythiaOptions->init();

    _description.reset(new po::options_description());
    _description->add(*_hiddenDescription)
        .add(*_muonInJetOptions->description())
        .add(*_pythiaOptions->description());
}

po::options_description *SolverInputOptions::description() const
{
    return _description.get();
}

void SolverInputOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    _muonInJetOptions->print(out);
    out << endl;

    _pythiaOptions->print(out);
    out << endl;

    out << "Solver Input Options" << endl;
    out << setw(25) << left << " [+] Input Type"
        << (_isData ? "data" : "monte-carlo") << endl;
}



void SolverInputOptions::optionDataIsSet(const bool &value)
{
    _isData = value;

    if (_delegate)
        _delegate->optionDataIsSet(value);
}
