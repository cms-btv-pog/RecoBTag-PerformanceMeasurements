/**
 * MultiStreamBuf
 * core 
 *
 * Created by Samvel Khalatian on Sep 10, 2010
 * Copyright 2010, All rights reserved
 */

#include "interface/MultiStreamBuffer.h"

using core::MultiStreamBuffer;

MultiStreamBuffer::MultiStreamBuffer()
{
}

void MultiStreamBuffer::add(std::streambuf *stream)
{
    _buffers.push_back(stream);
}

int MultiStreamBuffer::overflow(int c)
{
    for(BufferCollection::iterator buffer = _buffers.begin();
        _buffers.end() != buffer;
        ++buffer)
    {
        if (EOF == (*buffer)->sputc(c))
            return EOF;
    }

    return c;
}
