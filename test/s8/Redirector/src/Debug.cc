/**
 * Debug
 * core
 *
 * Created by Samvel Khalatian on Aug 30, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>

#include "interface/Debug.h"

using std::cerr;
using std::clog;
using std::cout;
using std::endl;

using core::Debug;

Debug::Debug():
    _cerrRedirector(cerr),
    _clogRedirector(clog),
    _coutRedirector(cout)
{
    // Supress Log by default
    //
    _clogRedirector.rdbuf(0);

    // Make COUT and CERR output into CLOG as well
    //
    _cerrMultiStream.reset(new MultiStreamBuffer());
    _cerrMultiStream->add(cerr.rdbuf());
    _cerrRedirector.rdbuf(_cerrMultiStream.get());

    _coutMultiStream.reset(new MultiStreamBuffer());
    _coutMultiStream->add(cout.rdbuf());
    _coutRedirector.rdbuf(_coutMultiStream.get());
}

void Debug::init(const std::string &filename)
try
{
    if (_debugFile.get())
    {
        cerr << "[Debug] debug file is already set. Can not redirect it to: "
            << filename << endl;

        return;
    }

    // Attemt to create output stream
    //
    _debugFile.reset(new std::ofstream(filename.c_str()));
    if (!_debugFile->is_open())
    {
        cerr << "[Debug] Failed to open output debug file: "
            << filename << endl;

        _debugFile.reset();
    }
    else
    {
        _clogRedirector.rdbuf(_debugFile->rdbuf());

        // Add debug file to the Multi Stream Buffers
        //
        _cerrMultiStream->add(_debugFile->rdbuf());
        _coutMultiStream->add(_debugFile->rdbuf());
    }
}
catch(...)
{
    // Failed to allocate new resource
    //
    cerr << "[Debug] Failed to allocate resources to open debuf file."
        << endl;
}
