/**
 * Converter
 * s8
 *
 * Created by Samvel Khalatian on Oct 1, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_CONVERTER
#define S8_CONVERTER

#include <string>

class TH2;

namespace s8
{
    class Event;
    class Muon;
    class Jet;
    class OperatingPoint;
    class Plots;

    // Generate S8Solver input histograms from S8Tuple
    //
    class Converter
    {
        public:
            Converter() throw();
            ~Converter() throw();

            bool run(const int &argc, char **argv);

        private:
            bool parseArguments(const int &argc, char **);

            void process();
            void analyze(const s8::Event *);

            struct Config
            {
                int         events;
                bool        isData;
                std::string tag;
                std::string output;
                std::string input;
            };

            Config         _config;

            Plots *_n;
            Plots *_p;
    };
}

#endif
