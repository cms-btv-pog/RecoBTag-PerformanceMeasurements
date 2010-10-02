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

            struct PlotGroup {
                PlotGroup(const std::string &, const std::string & = "");
                ~PlotGroup();

                void setOperatingPoint(const OperatingPoint &);

                void fill(const s8::Muon *, const s8::Jet *);

                void save() const;

                TH2 *all;
                TH2 *tag;

                private:
                    bool   _isInitializationFailed;
                    double _operatingPoint;
            };

            struct Plots
            {
                virtual ~Plots();

                virtual void setOperatingPoint(const OperatingPoint &) = 0;
                virtual void fill(const s8::Muon *, const s8::Jet *) = 0;
                virtual void save() const = 0;
            };

            struct NonFlavouredPlots : public Plots
            {
                NonFlavouredPlots(const std::string &);

                virtual void setOperatingPoint(const OperatingPoint &);
                virtual void fill(const s8::Muon *, const s8::Jet *);
                virtual void save() const;

                PlotGroup plots;
            };

            struct FlavouredPlots : public Plots
            {
                FlavouredPlots(const std::string &);

                virtual void setOperatingPoint(const OperatingPoint &);
                virtual void fill(const s8::Muon *, const s8::Jet *);
                virtual void save() const;

                PlotGroup b;
                PlotGroup cl;
            };

            Plots *_n;
            Plots *_p;
    };
}

#endif
