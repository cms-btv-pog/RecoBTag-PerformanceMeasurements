/**
 * MiscOptions
 * s8
 *
 * Created by Samvel Khalatian on Mar 9, 2011
 * Copyright 2011, All rights reserved
 */

#include <iomanip>
#include <stdexcept>
#include <string>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include "Option/interface/MiscOptionsDelegate.h"

#include "Option/interface/MiscOptions.h"

using std::runtime_error;
using std::string;

using s8::MiscOptions;
using s8::MiscOptionsDelegate;

MiscOptions::MiscOptions() throw()
{
    _delegate = 0;
}

MiscOptions::~MiscOptions() throw()
{
}

MiscOptionsDelegate *MiscOptions::delegate() const
{
    return _delegate;
}

void MiscOptions::setDelegate(MiscOptionsDelegate *delegate)
{
    _delegate = delegate;
}

void MiscOptions::init()
{
    _description.reset(new po::options_description("Misc Options"));
    _description->add_options()
        ("primary-vertices",
            po::value<string>()->notifier(
                boost::bind(&MiscOptions::optionPrimaryVerticesIsSet,
                            this, _1)),
            "Limit number of the primary vertices in event. Format: --primary-vertices=MIN..MAX")
    ;
}

po::options_description *MiscOptions::description() const
{
    return _description.get();
}

void MiscOptions::print(std::ostream &out) const
{
    using std::endl;
    using std::left;
    using std::setw;

    out << "Misc Options" << endl;
    out << setw(25) << left << " [+] Primary Vertices"
        << _primary_vertices << endl;
}



void MiscOptions::optionPrimaryVerticesIsSet(const string &value)
{
    if (value.empty())
        return;

    parse(_primary_vertices, value);

    if (_delegate)
        _delegate->optionPrimaryVerticesIsSet(_primary_vertices);
}
