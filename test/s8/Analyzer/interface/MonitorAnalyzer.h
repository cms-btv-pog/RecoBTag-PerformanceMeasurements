/**
 * MonitorAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MONITOR_ANALYZER
#define S8_MONITOR_ANALYZER

#include <memory>

#include <boost/shared_ptr.hpp>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/MonitorOptionsDelegate.h"
#include "Option/interface/PythiaOptionsDelegate.h"
#include "Utility/interface/Range.h"
#include "Utility/interface/MuonInJetDelegate.h"

class TFile;
class TH1;
class TH3;

namespace s8
{
    class MonitorDelta;
    class MonitorJet;
    class MonitorLepton;
    class MonitorMuonInJet;
    class MuonInJet;
    class S8Selector;
    class TaggerOperatingPoint;

    class MonitorAnalyzer: virtual public Analyzer,
                           virtual public MonitorOptionsDelegate,
                           virtual public MuonInJetDelegate
    {
        public:
            MonitorAnalyzer() throw();
            virtual ~MonitorAnalyzer() throw();

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
            virtual void optionSimulateTriggerIsSet(const bool &);
            virtual void optionReweightTriggerIsSet(const std::string &);

            // Misc options
            //
            virtual void optionPrimaryVerticesIsSet(const Range &);

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
            std::auto_ptr<TH1> _primary_vertices;

            std::auto_ptr<MonitorLepton> _monitorMuons;
            std::auto_ptr<MonitorJet> _monitorJets;
            std::auto_ptr<MonitorJet> _monitorLeadingJet;

            // (n) plots
            //
            std::auto_ptr<MonitorMuonInJet> _monitorNMuonInJet;
            std::auto_ptr<TH1> _n_primary_vertices;
            std::auto_ptr<TH1> _n_primary_vertices_z;
            bool _are_n_plots_filled;

            // (p) plots
            //
            std::auto_ptr<MonitorMuonInJet> _monitorPMuonInJet;

            Range _n_muon_pt;
            Range _n_jet_pt;
            Range _n_jet_eta;

            S8Selector    *_s8_selector;

            boost::shared_ptr<TFile> _reweight_file;
            TH3 *_reweights;

            const Event *_event;
    };
}

#endif
