/**
 * Monitor Trigger analyzer
 * 
 *
 * Created by Samvel Khalatian on Feb 25, 2011
 * Copyright 2011, All rights reserved
 */

#include <iostream>
#include <memory>
#include <stdexcept>

#include "Controller/interface/MonitorTriggerController.h"

using std::bad_alloc;
using std::cerr;
using std::endl;

int main(int argc, char **argv) 
try
{
    std::auto_ptr<core::AppController> app(new s8::MonitorTriggerController());
    if (app->init(argc, argv))
        return app->run() ? 0 : 1;

    return 0;
}
catch(const bad_alloc &error)
{
    cerr << "[error] Failed to instantiate Controller. Insufficient memory?" << endl;
    cerr << endl;
    cerr << "system error message:" << endl;
    cerr << endl;
    cerr << error.what() << endl;

    return 1;
}
catch(...)
{
    cerr << "[error] unknown exception was raised." << endl;

    return 1;
}
