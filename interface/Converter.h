/**
 * Converter
 * 
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

    // Generate S8Solver input histograms from S8Tuple
    //
    class Converter
    {
        public:
            Converter() throw();
            ~Converter() throw();

            bool run(const int argc, char *argv[]);

        private:
            void process(const std::string &file);
            void analyze(const s8::Event *);

            struct Plots {
                Plots(const std::string &, const std::string &);
                ~Plots();

                void setOperatingPoint(const OperatingPoint &);

                void fill(const s8::Muon *, const s8::Jet *);

                void save() const;

                TH2 *all;
                TH2 *tag;

                private:
                    bool   _isInitializationFailed;
                    double _operatingPoint;
            };

            struct FlavouredPlots {
                FlavouredPlots(const std::string &);

                void setOperatingPoint(const OperatingPoint &);

                void fill(const s8::Muon *, const s8::Jet *);

                void save() const;

                Plots b;
                Plots cl;
            };

            FlavouredPlots _n;
            FlavouredPlots _p;
    };
}

#endif
