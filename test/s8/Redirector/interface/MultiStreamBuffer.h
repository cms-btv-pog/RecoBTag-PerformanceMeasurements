/**
 * MultiStreamBuffer
 * core 
 *
 * Created by Samvel Khalatian on Sep 10, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_MULTI_STREAM_BUFFER
#define CORE_MULTI_STREAM_BUFFER

#include <vector>
#include <streambuf>

namespace core
{
    class MultiStreamBuffer : public std::streambuf
    {
        public:
            MultiStreamBuffer();

            void add(std::streambuf *);

        protected:
            virtual int overflow(int c);

        private:
            typedef std::vector<std::streambuf *> BufferCollection;

            BufferCollection _buffers;
    };
}

#endif
