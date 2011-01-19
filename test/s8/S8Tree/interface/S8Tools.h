/**
 * Tools
 * s8 
 *
 * Created by Samvel Khalatian on Nov 9, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TOOLS
#define S8_TOOLS

#include <string>

#include <boost/functional/hash.hpp>

namespace s8
{
    namespace tools
    {
        typedef unsigned int Hash;

        typedef boost::hash<std::string> HashGenerator;

        HashGenerator make_hash;
    }
}

#endif
