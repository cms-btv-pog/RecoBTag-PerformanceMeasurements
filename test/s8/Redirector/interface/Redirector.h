/**
 * Redirector
 * core 
 *
 * Created by Samvel Khalatian on Sep 10, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_REDIRECTOR
#define CORE_REDIRECTOR

#include <iostream>
#include <streambuf>

namespace core
{
    class Redirector
    {
        public:
            Redirector(std::ostream &);
            ~Redirector();

            std::streambuf *rdbuf(std::streambuf *);

        private:
            std::ostream   &_stream;
            std::streambuf *_buffer;
    };
}

#endif
