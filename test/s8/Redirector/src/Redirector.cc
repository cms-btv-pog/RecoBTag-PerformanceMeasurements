/**
 * Redirector
 * core 
 *
 * Created by Samvel Khalatian on Sep 10, 2010
 * Copyright 2010, All rights reserved
 */

#include "interface/Redirector.h"

using core::Redirector;

Redirector::Redirector(std::ostream &stream):
    _stream(stream)
{
    _buffer = 0;
}

Redirector::~Redirector()
{
    if (!_buffer)
        return;

    _stream.rdbuf(_buffer);
}

std::streambuf *Redirector::rdbuf(std::streambuf *buffer)
{
    if (!_buffer)
    {
        _buffer = _stream.rdbuf();
    }

    return _stream.rdbuf(buffer);
}
