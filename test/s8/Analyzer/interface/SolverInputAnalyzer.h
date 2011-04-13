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

#include <boost/shared_ptr.hpp>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/SolverInputOptionsDelegate.h"
#include "Utility/interface/Range.h"
#include "Utility/interface/MuonInJetDelegate.h"

class TFile;
class TH3;

namespace s8
{
    class LeptonInJetDiscriminator;
    class MuonInJet;
    class S8Selector;
    class SolverLeptonInJetPlots;
    class TaggerOperatingPoint;
    class Trigger;

    class SolverInputAnalyzer: virtual public Analyzer,
                               virtual public SolverInputOptionsDelegate,
                               virtual public MuonInJetDelegate
    {
        public:
            SolverInputAnalyzer() throw();
            virtual ~SolverInputAnalyzer() throw();

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

            // MuonInJetDelegate interface
            //
            virtual bool shouldSkipMuonInJetPlusAwayJet(const Lepton *,
                                                        const Jet *);

            virtual void muonIsInJetPlusAwayJet(const Lepton *,
                                                const Jet *);

            virtual void muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                      const Jet *);

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
            std::auto_ptr<TaggerOperatingPoint> _muonInJetTagger;
            std::auto_ptr<TaggerOperatingPoint> _awayJetTagger;

            std::auto_ptr<MuonInJet> _muonInJet;

            Range _n_muon_pt;
            Range _n_jet_pt;
            Range _n_jet_eta;

            S8Selector    *_s8_selector;

            std::auto_ptr<LeptonInJetDiscriminator> _leptonInJetDiscriminator;
            std::auto_ptr<SolverLeptonInJetPlots>   _solverNPlots;
            std::auto_ptr<SolverLeptonInJetPlots>   _solverPPlots;

            boost::shared_ptr<TFile> _reweight_file;
            TH3 *_reweights;
            int _n_primary_vertices;
    };
}

#endif
