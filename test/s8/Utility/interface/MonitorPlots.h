/**
 * MonitorPlots
 * s8
 *
 * Created by Samvel Khalatian on Oct 15, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MONITOR_PLOTS
#define S8_MONITOR_PLOTS

#include <memory>
#include <string>

class TDirectory;
class TH1;
class TH3;
class TLorentzVector;

namespace s8
{
    class Jet;
    class Lepton;

    class Monitor
    {
        public:
            Monitor();
            virtual ~Monitor() throw();

            virtual void save(TDirectory *) = 0;
    };

    class MonitorBase: public Monitor
    {
        public:
            enum Plot { PT, ETA, PHI };

            MonitorBase(const std::string &);
            virtual ~MonitorBase() throw();

            virtual void save(TDirectory *);

            TH1 *plot(const Plot &);

        protected:
            void fill(const double &weight, const Plot &, const double &);

        private:
            std::auto_ptr<TH1> _pt;
            std::auto_ptr<TH1> _eta;
            std::auto_ptr<TH1> _phi;
    };

    class MonitorDelta: public Monitor
    {
        public:
            enum Plot { PTREL, DELTA_R, DELTA_PHI, DELTA_ETA };

            MonitorDelta(const std::string &);
            virtual ~MonitorDelta() throw();

            virtual void save(TDirectory *);

            void fill(const double &weight, const TLorentzVector *, const TLorentzVector *);

            TH1 *plot(const Plot &);

        private:
            std::auto_ptr<TH1> _ptrel;
            std::auto_ptr<TH1> _deltaR;
            std::auto_ptr<TH1> _deltaPhi;
            std::auto_ptr<TH1> _deltaEta;
    };

    class MonitorLepton: public MonitorBase
    {
        public:
            MonitorLepton(const std::string &prefix);

            void fill(const double &weight, const Lepton *);
    };

    class MonitorJet: public MonitorBase
    {
        public:
            MonitorJet(const std::string &prefix);
            virtual ~MonitorJet() throw();

            void fill(const double &weight, const Jet *);
            void fillDiscriminator(const double &weight, const double &value);

            virtual void save(TDirectory *);

        private:
            std::auto_ptr<TH1> _discriminator;
    };

    class MonitorMuonInJet: public Monitor
    {
        public:
            enum Plot { MUON, JET, DELTA };

            // Use: 'n' or 'p' as prefix
            //
            MonitorMuonInJet(const std::string &prefix);
            virtual ~MonitorMuonInJet() throw();

            MonitorJet *jet();
            MonitorDelta *delta();

            void fill(const double &weight, const Lepton *, const Jet *,
                    const int &npv = 0);

            virtual void save(TDirectory *);

        private:
            MonitorLepton _muon;
            MonitorJet    _jet;
            MonitorDelta  _delta;
            std::auto_ptr<TH3> _pt_eta_pv;
    };
}

#endif
