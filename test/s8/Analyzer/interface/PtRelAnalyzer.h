/**
 * PtRelAnalyzer
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PTREL_ANALYZER
#define S8_PTREL_ANALYZER

#include <memory>

#include "Analyzer/interface/Analyzer.h"
#include "Option/interface/SolverInputOptionsDelegate.h"
#include "Utility/interface/MuonInJetDelegate.h"

class TH2;

namespace s8
{
    class Jet;
    class Lepton;
    class MuonInJet;
    class TaggerOperatingPoint;

    class PtRelAnalyzer: public Analyzer,
                         public SolverInputOptionsDelegate,
                         public MuonInJetDelegate
    {
        public:
            PtRelAnalyzer() throw();
            virtual ~PtRelAnalyzer() throw();

            // Analyzer interface
            //
            virtual void init();

            virtual void eventDidLoad(const Event *);
            virtual void print(std::ostream &) const;
            virtual void save(TDirectory *) const;

            // SolverInputOptionsDelegate interface
            //
            virtual void optionAwayTagIsSet(const std::string &);
            virtual void optionDataIsSet(const bool &);

            // MuonInJetDelegate interface
            //
            virtual void muonIsInJetPlusAwayJet(const Lepton *, const Jet *);

        private:
            std::auto_ptr<TaggerOperatingPoint> _awayJetTagger;

            std::auto_ptr<MuonInJet> _muonInJet;

            std::auto_ptr<TH2> _ptrel;
            std::auto_ptr<TH2> _ptrel_b;
            std::auto_ptr<TH2> _ptrel_c;
            std::auto_ptr<TH2> _ptrel_cl;
            std::auto_ptr<TH2> _ptrel_l;
            std::auto_ptr<TH2> _ptrel_g;

            bool _isData;
    };
}

#endif
