/**
 * SolverInputAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_SOLVER_INPUT_ANALYZER
#define S8_SOLVER_INPUT_ANALYZER

#include <memory>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/SolverInputOptionsDelegate.h"
#include "Utility/interface/MuonInJetDelegate.h"

namespace s8
{
    class Jet;
    class Lepton;
    class LeptonInJetDiscriminator;
    class MuonInJet;
    class SolverLeptonInJetPlots;
    class TaggerOperatingPoint;

    class SolverInputAnalyzer: public Analyzer,
                               public SolverInputOptionsDelegate,
                               public MuonInJetDelegate
    {
        public:
            SolverInputAnalyzer() throw();
            virtual ~SolverInputAnalyzer() throw();

            // Analyzer interface
            //
            virtual void init();

            virtual void eventDidLoad(const Event *);

            virtual void save(TDirectory *) const;

            // SolverInputOptionsDelegate interface
            //
            virtual void optionDataIsSet(const bool &);
            virtual void optionGluonSplittingIsSet(const GluonSplitting &);
            virtual void optionPtHatIsSet(const int &, const int &);

            // MuonInJetOptionsDelegate interface
            //
            virtual void optionTagIsSet(const std::string &);
            virtual void optionAwayTagIsSet(const std::string &);
            virtual void optionMuonPtIsSet(const double &);

            // MuonInJetDelegate interface
            //
            virtual void muonIsInJetPlusAwayJet(const Lepton *,
                                                const Jet *);

            virtual void muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                      const Jet *);

        private:
            std::auto_ptr<TaggerOperatingPoint> _muonInJetTagger;
            std::auto_ptr<TaggerOperatingPoint> _awayJetTagger;

            std::auto_ptr<MuonInJet> _muonInJet;

            GluonSplitting _gluonSplitting;
            int            _minPtHat;
            int            _maxPtHat;
            bool           _otherEvent; // Count every second event for
                                        // gluon splitting

            std::auto_ptr<LeptonInJetDiscriminator> _leptonInJetDiscriminator;
            std::auto_ptr<SolverLeptonInJetPlots>   _solverNPlots;
            std::auto_ptr<SolverLeptonInJetPlots>   _solverPPlots;
    };
}

#endif
