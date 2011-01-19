/**
 * SolverInputJetPtAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_SOLVER_INPUT_JETPT_ANALYZER
#define S8_SOLVER_INPUT_JETPT_ANALYZER

#include <memory>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/SolverInputOptionsDelegate.h"
#include "Utility/interface/TwoMuonsInJetDelegate.h"

namespace s8
{
    class Jet;
    class Lepton;
    class LeptonInJetDiscriminator;
    class TwoMuonsInJet;
    class SolverLeptonInJetPlots;
    class TaggerOperatingPoint;

    class SolverInputJetPtAnalyzer: public Analyzer,
                                    public SolverInputOptionsDelegate,
                                    public TwoMuonsInJetDelegate
    {
        public:
            SolverInputJetPtAnalyzer() throw();
            virtual ~SolverInputJetPtAnalyzer() throw();

            // Analyzer interface
            //
            virtual void init();

            virtual void eventDidLoad(const Event *);

            virtual void save(TDirectory *) const;

            // SolverInputOptionsDelegate interface
            //
            virtual void optionDataIsSet(const bool &);

            // MuonInJetOptionsDelegate interface
            //
            virtual void optionTagIsSet(const std::string &);
            virtual void optionAwayTagIsSet(const std::string &);

            // TwoMuonsInJetDelegate interface
            //
            virtual void muonIsInJetPlusAwayJet(const Lepton *,
                                                const Lepton *,
                                                const Jet *);

            virtual void muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                      const Lepton *,
                                                      const Jet *);

        private:
            std::auto_ptr<TaggerOperatingPoint> _muonInJetTagger;
            std::auto_ptr<TaggerOperatingPoint> _awayJetTagger;

            std::auto_ptr<TwoMuonsInJet> _muonInJet;

            std::auto_ptr<LeptonInJetDiscriminator> _leptonInJetDiscriminator;
            std::auto_ptr<SolverLeptonInJetPlots>   _solverNPlots;
            std::auto_ptr<SolverLeptonInJetPlots>   _solverPPlots;
    };
}

#endif
