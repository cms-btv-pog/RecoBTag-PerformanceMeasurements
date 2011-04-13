/**
 * MonitorTriggerAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Feb 25, 2011
 * Copyright 2011, All rights reserved
 */

#ifndef S8_MONITOR_TRIGGER_ANALYZER
#define S8_MONITOR_TRIGGER_ANALYZER

#include <memory>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/MonitorOptionsDelegate.h"
#include "Option/interface/PythiaOptionsDelegate.h"
#include "S8Tree/interface/S8Trigger.h"
#include "Selector/interface/TriggerNamePredicator.h"
#include "Utility/interface/Range.h"
#include "Utility/interface/MuonInJetDelegate.h"

class TH1;

namespace s8
{
    class MonitorDelta;
    class MonitorJet;
    class MonitorLepton;
    class MuonInJet;
    class S8Selector;
    class TaggerOperatingPoint;

    class MonitorTriggerAnalyzer: public Analyzer,
                                  public MonitorOptionsDelegate,
                                  public MuonInJetDelegate
    {
        public:
            MonitorTriggerAnalyzer() throw();
            virtual ~MonitorTriggerAnalyzer() throw();

            // Analyzer interface
            //
            virtual void init();

            virtual void treeDidLoad(const TreeInfo *, const TriggerCenter *);
            virtual void eventDidLoad(const Event *);
            virtual void print(std::ostream &) const;
            virtual void save(TDirectory *) const;

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
            virtual void optionUseTriggerPrescaleIsSet(const bool &);

        private:
            void saveGenericPlots(TDirectory *) const;
            void saveNPlots(TDirectory *) const;
            void savePPlots(TDirectory *) const;

            std::auto_ptr<MuonInJet> _muonInJet;

            std::auto_ptr<TaggerOperatingPoint> _muonInJetTagger;
            std::auto_ptr<TaggerOperatingPoint> _awayJetTagger;

            // Generic plots
            //
            std::auto_ptr<TH1> _muons;
            std::auto_ptr<TH1> _jets;
            std::auto_ptr<TH1> _pthat;

            std::auto_ptr<MonitorLepton> _monitorMuons;
            std::auto_ptr<MonitorJet> _monitorJets;
            std::auto_ptr<MonitorJet> _monitorLeadingJet;

            // (n) plots
            //
            std::auto_ptr<MonitorLepton> _monitorNMuons;
            std::auto_ptr<MonitorJet>    _monitorNJets;
            std::auto_ptr<MonitorDelta>  _monitorNDelta;

            // (p) plots
            //
            std::auto_ptr<MonitorLepton> _monitorPMuons;
            std::auto_ptr<MonitorJet>    _monitorPJets;
            std::auto_ptr<MonitorDelta>  _monitorPDelta;

            Range _n_muon_pt;
            Range _n_jet_pt;
            Range _n_jet_eta;

            S8Selector    *_s8_selector;

            bool _use_trigger_prescale;
            int _prescale;
            Trigger _trigger;
            TriggerNamePredicator _trigger_predicator;
    };
}

#endif
