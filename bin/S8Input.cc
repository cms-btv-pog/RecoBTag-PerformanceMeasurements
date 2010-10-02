/**
 * Converter
 * 
 *
 * Created by Samvel Khalatian on Oct 1, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <memory>
#include <stdexcept>

#include "RecoBTag/PerformanceMeasurements/interface/Converter.h"

using std::auto_ptr;
using std::cerr;
using std::exception;
using std::endl;

using s8::Converter;

int main(int argc, char *argv[])
try
{
    auto_ptr<Converter> converter(new Converter());

    return converter->run(argc,argv) ?
        0 :
        1;
}
catch(const exception &error)
{
    cerr << "[error] " << error.what() << endl;

    return 0;
}
catch(...)
{
    cerr << "unknown error. Abnormal exit." << endl;

    return 0;
}
