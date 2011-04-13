/**
 * SolverLeptonInJetPlots
 * s8 
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_SOLVER_MUON_IN_JET_PLOTS
#define S8_SOLVER_MUON_IN_JET_PLOTS

#include <memory>
#include <string>

class TDirectory;
class TH2;

namespace s8
{
    class Jet;
    class Lepton;
    class TaggerOperatingPoint;
    class LeptonInJetDiscriminator;

    // plot namespace is designed solely for internal purposes. User should be
    // working with s8::SolverLeptonInJetPlots class
    //
    namespace plot
    {
        // Generic interface for all plots
        //
        class LeptonInJetPlots
        {
            public:
                virtual ~LeptonInJetPlots();

                // Plots should be booked in init
                //
                virtual void init() = 0;

                virtual void setDiscriminator(LeptonInJetDiscriminator *) = 0;

                // Group is responsible for filling plots
                //
                virtual void fill(const Lepton *, const Jet *,
                                  const double &weight = 1) = 0;

                // Each group knows how to save itself
                //
                virtual void save(TDirectory *) const = 0;
        };

        // Generic interface for all groups
        //
        class LeptonInJetGroup
        {
            public:
                virtual ~LeptonInJetGroup();

                // TaggerOperatingPoint is used for b-Discriminant
                //
                virtual void setTaggerOperatingPoint(const TaggerOperatingPoint *) = 0;
        };

        // Group of plots that has 2D Pt, Eta, etc. vs PtRel distributions,e.g.
        //
        //      PtRel vs Pt
        //      PtRel vs Eta
        //      PtRel vs etc.
        //
        //  where:
        //
        //  PtRel is the length of the perpendicular vector of the muon momentum
        //        with respect to (muon+jet) momentum
        //
        //  Pt    is the corresponding Jet Pt
        //
        //  Eta   is the corresponding Jet Eta
        //
        class LeptonInJetPlotCategory: public LeptonInJetPlots
        {
            public:
                // Prefix and Suffix
                //
                LeptonInJetPlotCategory(const std::string &, const std::string &);
                virtual ~LeptonInJetPlotCategory();

                virtual void init();

                virtual void setDiscriminator(LeptonInJetDiscriminator *);

                virtual void fill(const Lepton *, const Jet *,
                                  const double &weight = 1);
                virtual void save(TDirectory *) const;

            private:
                bool _didNotInitialize;

                std::string _prefix;
                std::string _suffix;

                LeptonInJetDiscriminator *_discriminator;

                std::auto_ptr<TH2>          _pt;
                std::auto_ptr<TH2>          _eta;
                std::auto_ptr<TH2>          _phi;
        };

        // All and b-tagged jet categories.
        //
        class LeptonInJetPlotGroup: public LeptonInJetPlots,
                                    public LeptonInJetGroup
        {
            public:
                LeptonInJetPlotGroup(const std::string &, const std::string & = "");
                virtual ~LeptonInJetPlotGroup();

                virtual void init();

                virtual void setDiscriminator(LeptonInJetDiscriminator *);

                virtual void setTaggerOperatingPoint(const TaggerOperatingPoint *);
                virtual void fill(const Lepton *, const Jet *,
                                  const double &weight = 1);
                virtual void save(TDirectory *) const;

            private:
                const TaggerOperatingPoint *_operatingPoint;

                std::auto_ptr<LeptonInJetPlotCategory> _all;
                std::auto_ptr<LeptonInJetPlotCategory> _tag;
        };

        // One PlotGroup for all jet flavours. Given set of plots should be
        // used by Data input.
        //
        class NonFlavouredLeptonInJetPlotGroup: public LeptonInJetPlots,
                                                public LeptonInJetGroup
        {
            public:
                NonFlavouredLeptonInJetPlotGroup(const std::string &);

                virtual void init();

                virtual void setDiscriminator(LeptonInJetDiscriminator *);

                virtual void setTaggerOperatingPoint(const TaggerOperatingPoint *);
                virtual void fill(const Lepton *, const Jet *,
                                  const double &weight = 1);
                virtual void save(TDirectory *) const;

            private:
                std::auto_ptr<LeptonInJetPlotGroup> _plots;
        };

        // Separate PlotGroup for each jet flavour
        //
        class FlavouredLeptonInJetPlotGroup: public LeptonInJetPlots,
                                             public LeptonInJetGroup
        {
            public:
                FlavouredLeptonInJetPlotGroup(const std::string &);

                virtual void init();

                virtual void setDiscriminator(LeptonInJetDiscriminator *);

                virtual void setTaggerOperatingPoint(const TaggerOperatingPoint *);
                virtual void fill(const Lepton *, const Jet *,
                                  const double &weight = 1);
                virtual void save(TDirectory *) const;

            private:
                std::auto_ptr<LeptonInJetPlotGroup> _b;
                std::auto_ptr<LeptonInJetPlotGroup> _cl;
        };
    }

    // Combination of flavoured and non-flavoured plots. This set of plots
    // should be used by Monte-Carlo inputs.
    //
    class SolverLeptonInJetPlots: public plot::LeptonInJetPlots,
                                  public plot::LeptonInJetGroup
    {
        public:
            SolverLeptonInJetPlots(const std::string &);

            virtual void init();

            bool isFlavoured() const;
            void setIsFlavoured(const bool &);

            virtual void setDiscriminator(LeptonInJetDiscriminator *);

            virtual void setTaggerOperatingPoint(const TaggerOperatingPoint *);
            virtual void fill(const Lepton *, const Jet *,
                              const double &weight = 1);
            virtual void save(TDirectory *) const;

        private:
            bool _isFlavoured;

            std::auto_ptr<plot::FlavouredLeptonInJetPlotGroup>    _flavoured;
            std::auto_ptr<plot::NonFlavouredLeptonInJetPlotGroup> _nonFlavoured;
    };
}

#endif
