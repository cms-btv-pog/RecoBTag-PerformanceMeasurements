/**
 * Debug
 * core 
 *
 * Created by Samvel Khalatian on Aug 30, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_DEBUG
#define CORE_DEBUG

#include <fstream>
#include <memory>

#include "interface/MultiStreamBuffer.h"
#include "interface/Redirector.h"

namespace core
{
    class Debug
    {
        public:
            Debug();

            void init(const std::string &filename);

        private:
            // WARNING: Do not change order ! Destruction order is important:
            //  1. Close file if open
            //  2. Restore Output streams
            //  3. Destroy Multi Stream Buffers
            //
            std::auto_ptr<std::ofstream> _debugFile;

            std::auto_ptr<MultiStreamBuffer> _cerrMultiStream;
            std::auto_ptr<MultiStreamBuffer> _coutMultiStream;

            Redirector _cerrRedirector;
            Redirector _clogRedirector;
            Redirector _coutRedirector;
    };
}

#endif
