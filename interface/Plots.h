/**
 * Plots
 * s8 
 *
 * Created by Samvel Khalatian on Oct 3, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PLOTS
#define S8_PLOTS

#include <string>

class TDirectory;
class TH2;

namespace s8
{
    class Jet;
    class Muon;
    class OperatingPoint;

    struct Plots
    {
        virtual ~Plots();

        virtual void setOperatingPoint(const OperatingPoint &) = 0;
        virtual void fill(const Muon *, const Jet *) = 0;
        virtual void save(TDirectory *) const = 0;
    };

    class PlotGroup : public Plots
    {
        public:
            PlotGroup(const std::string &, const std::string & = "");
            ~PlotGroup();

            void setOperatingPoint(const OperatingPoint &);
            void fill(const Muon *, const Jet *);
            void save(TDirectory *) const;

        private:
            bool   _isInitializationFailed;
            double _operatingPoint;

            TH2 *_allPt;
            TH2 *_tagPt;

            TH2 *_allEta;
            TH2 *_tagEta;
    };

    class NonFlavouredPlots : public Plots
    {
        public:
            NonFlavouredPlots(const std::string &);

            virtual void setOperatingPoint(const OperatingPoint &);
            virtual void fill(const Muon *, const Jet *);
            virtual void save(TDirectory *) const;

        private:
            PlotGroup _plots;
    };

    class FlavouredPlots : public Plots
    {
        public:
            FlavouredPlots(const std::string &);

            virtual void setOperatingPoint(const OperatingPoint &);
            virtual void fill(const Muon *, const Jet *);
            virtual void save(TDirectory *) const;

        private:
            PlotGroup _b;
            PlotGroup _cl;
    };

    class CombinedPlots : public Plots
    {
        public:
            CombinedPlots(const std::string &);

            virtual void setOperatingPoint(const OperatingPoint &);
            virtual void fill(const Muon *, const Jet *);
            virtual void save(TDirectory *) const;

        private:
            FlavouredPlots    _flavoured;
            NonFlavouredPlots _nonFlavoured;
    };
}

#endif
