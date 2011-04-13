/**
 * MonitorAnalyzer + SolverInput
 * s8
 *
 * Created by Samvel Khalatian on Feb 14, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MONITOR_SOLVER_INPUT_ANALYZER
#define S8_MONITOR_SOLVER_INPUT_ANALYZER

#include <boost/shared_ptr.hpp>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/MonitorOptionsDelegate.h"
#include "Option/interface/SolverInputOptionsDelegate.h"

namespace s8
{
    class MonitorAnalyzer;
    class SolverInputAnalyzer;

    class MonitorSolverInputAnalyzer: public Analyzer,
                                      public MonitorOptionsDelegate,
                                      public SolverInputOptionsDelegate
    {
        public:
            MonitorSolverInputAnalyzer() throw();
            virtual ~MonitorSolverInputAnalyzer() throw();

            // Analyzer interface
            //
            virtual void init();

            virtual void treeDidLoad(const TreeInfo *, const TriggerCenter *);
            virtual void eventDidLoad(const Event *);
            virtual void print(std::ostream &) const;
            virtual void save(TDirectory *) const;

            // SolverInputOptionsDelegate interface
            //
            virtual void optionDataIsSet(const bool &);

            // MuonInJet options
            //
            virtual void optionTagIsSet(const std::string &);
            virtual void optionAwayTagIsSet(const std::string &);
            virtual void optionMuonPtIsSet(const Range &);
            virtual void optionJetPtIsSet(const Range &);
            virtual void optionJetEtaIsSet(const Range &);

            // PythiaOptionsDelegate interface
            //
            virtual void optionGluonSplittingIsSet(const GluonSplitting &);
            virtual void optionPtHatIsSet(const Range &);

            // Trigger options
            //
            virtual void optionTriggerIsSet(const Trigger &);
            virtual void optionSimulateTriggerIsSet(const bool &);
            virtual void optionReweightTriggerIsSet(const std::string &);

            // Misc options
            //
            virtual void optionPrimaryVerticesIsSet(const Range &);

        private:
            boost::shared_ptr<MonitorAnalyzer> _monitor;
            boost::shared_ptr<SolverInputAnalyzer> _solver;
    };
}

#endif
