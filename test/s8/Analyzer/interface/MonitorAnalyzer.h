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

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/MonitorOptionsDelegate.h"
#include "Option/interface/PythiaOptionsDelegate.h"
#include "Utility/interface/MuonInJetDelegate.h"

class TH1;

namespace s8
{
    class MonitorDelta;
    class MonitorJet;
    class MonitorLepton;
    class MuonInJet;
    class TaggerOperatingPoint;

    class MonitorAnalyzer: public Analyzer,
                           public MonitorOptionsDelegate,
                           public MuonInJetDelegate
    {
        public:
            MonitorAnalyzer() throw();
            virtual ~MonitorAnalyzer() throw();

            // Analyzer interface
            //
            virtual void init();

            virtual void eventDidLoad(const Event *);

            virtual void save(TDirectory *) const;

            // MuonInJetOptionsDelegate interface
            //
            virtual void optionTagIsSet(const std::string &);
            virtual void optionAwayTagIsSet(const std::string &);

            // MuonInJetDelegate interface
            //
            virtual void muonIsInJetPlusAwayJet(const Lepton *,
                                                const Jet *);

            virtual void muonIsInJetPlusTaggedAwayJet(const Lepton *,
                                                      const Jet *);

            // PythiaOptionsDelegate interface
            //
            virtual void optionGluonSplittingIsSet(const GluonSplitting &);
            virtual void optionPtHatIsSet(const int &min, const int &max);

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

            int _minPtHat;
            int _maxPtHat;
    };
}

#endif
