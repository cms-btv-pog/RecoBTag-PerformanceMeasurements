#include <iostream>
#include <string>
#include <utility>

#include <boost/functional/hash.hpp>

#include <TClass.h>

#include "interface/S8TreeInfo.h"

using std::cout;
using std::endl;
using std::make_pair;
using std::string;

using s8::Version;

int main(int argc, char **argv)
{
    Version version1;
    version1.set(make_pair(10, 2));

    cout << "Version1: " << version1 << endl;

    Version version2;
    version2.set(make_pair(10, 11));

    cout << "Version2: " << version2 << endl;

    cout << "Version1 "
        << (version1 == version2 ? "==" : "!=") << " Version2" << endl;

    cout << "Version1 "
        << (version1 < version2 ? "<" : ">") << " Version2" << endl;

    cout << "Version1 "
        << (version1 > version2 ? ">" : "<") << " Version2" << endl;

    cout << "Version1 "
        << (version1 != version2 ? "!=" : "==") << " Version2" << endl;

    return 1;
}
