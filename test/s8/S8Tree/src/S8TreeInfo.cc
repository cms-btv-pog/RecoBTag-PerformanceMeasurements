/**
 * TreeInfo
 * s8 
 *
 * Created by Samvel Khalatian on Oct 22, 2010
 * Copyright 2010, All rights reserved
 */

#include <sstream>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include "interface/S8TreeInfo.h"

using std::make_pair;
using std::ostream;
using std::runtime_error;

using boost::lexical_cast;

using s8::TreeInfo;
using s8::Version;

Version::Version() throw():
    _version(make_pair(0, 0))
{
}

Version::VersionPair Version::operator()() const
{
    return _version;
}

void Version::set(const VersionPair &version)
{
    if (0 > version.first ||
        0 > version.second)

        throw runtime_error("Failed to set Version: version can not be negative");

    _version = version;
}

std::ostream &s8::operator<<(std::ostream &out, const Version &version)
{
    return out << version().first << "." << version().second;
}

bool s8::operator==(const Version &v1, const Version &v2)
{
    return v1() == v2();
}

bool s8::operator!=(const Version &v1, const Version &v2)
{
    return v1() != v2();
}

bool s8::operator <(const Version &v1, const Version &v2)
{
    return v1() < v2();
}

bool s8::operator<=(const Version &v1, const Version &v2)
{
    return v1() <= v2();
}

bool s8::operator >(const Version &v1, const Version &v2)
{
    return v1() > v2();
}

bool s8::operator>=(const Version &v1, const Version &v2)
{
    return v1() >= v2();
}



TreeInfo::TreeInfo() throw()
{
    _version.set(make_pair(5, 0));
}

Version TreeInfo::version() const
{
    return _version;
}

void TreeInfo::merge(const TreeInfo &treeInfo)
{
    if (treeInfo.version() != version())
    {
        std::ostringstream message("Can not merge Trees as they are of different versions: ");
        message << 'v' << version() << " < v" << treeInfo.version();

        throw runtime_error(message.str());
    }
}

ClassImp(TreeInfo)
