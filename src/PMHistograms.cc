#include <map>
#include <vector>
#include <string>

#include "RecoBTag/PerformanceMeasurements/interface/PMHistograms.h"

#include "TH2F.h"
#include "TString.h"
#include "TH3F.h"

using std::string;

//______________
void PMHistograms::Add()
{
    std::map<std::string, int>  quark_color;
    quark_color[""] = 1;
    quark_color["b"] = 2;
    quark_color["c"] = 3;
    quark_color["uds"] = 4;
    quark_color["g"] = 6;

    TString           ftagger;
    TString           flevel;
    TString           fAwaytagger;
    TString           fAwaylevel;
    std::map< TString, float > fTrackCountingMap;
    TAxis             fJetPtAxis;
    TAxis             fJetEtaAxis;
    TAxis             fCorrPtAxis;
    TAxis             fCorrEtaAxis;

    ftagger = "TrackCounting";
    flevel  = "Loose";
    fAwaytagger = "TrackCounting";
    fAwaylevel = "Loose";

    fTrackCountingMap["Loose"]  = 1.7; // use TC2:high eff.
    fTrackCountingMap["Medium"] = 3.3; // use TC2:high eff.
    fTrackCountingMap["Tight"]  = 3.4;

	// const int nptarray = 7;
	const int nptarray = 4;
    const int netaarray = 5;
    const int nphiarray = 8;
	const int nptrelarray = 51;
    //const int ncorrptarray = 3;
    //const int ncorretaarray = 5;
//    Double_t jetptbins[nptarray] = {10.,20.,30.,50., 70, 100., 230.};
    Double_t jetptbins[nptarray] = {30.,50.,80.,230.};
    Double_t jetetabins[netaarray] = {0.0,0.5,1.0,1.5,2.5};
    Double_t jetphibins[nphiarray] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5};
	Double_t ptrelbins[nptrelarray] = {0.,0.1,0.2,0.3,0.4,
									   0.5,0.6,0.7,0.8,0.9,
									   1.0,1.1,1.2,1.3,1.4,
									   1.5,1.6,1.7,1.8,1.9,
									   2.0,2.1,2.2,2.3,2.4,
									   2.5,2.6,2.7,2.8,2.9,
									   3.0,3.1,3.2,3.3,3.4,
									   3.5,3.6,3.7,3.8,3.9,
									   4.0,4.1,4.2,4.3,4.4,
									   4.5,4.6,4.7,4.8,4.9,
									   5.0};
	
    //Double_t corrptbins[ncorrptarray] = {20.,40.,60.,80.,230.};
    //Double_t corretabins[ncorrptarray] = {0.,0.5,1.,1.5,2.5};

    int nptbins = nptarray -1;
    //const Double_t *jetptbins = (fJetPtAxis.GetXbins())->GetArray();
    int netabins = netaarray -1;
    int nphibins = nphiarray - 1;
    int nptrelbins = nptrelarray -1;
    //const Double_t *jetetabins = (fJetEtaAxis.GetXbins())->GetArray();
    // int ncorrptbins = fCorrPtAxis.GetNbins();
    // const Double_t *corrptbins = (fCorrPtAxis.GetXbins())->GetArray();
    // int ncorretabins = fCorrEtaAxis.GetNbins();
    // const Double_t *corretabins = (fCorrEtaAxis.GetXbins())->GetArray();


	// Pt
	fstore->add( new TH3F("n3_pT","MuTag pT,eta vs pTrel",nptbins,jetptbins,nptrelbins,ptrelbins,netabins,jetetabins), "muon_in_jet" );

    using std::vector;
    using std::map;

    typedef vector<string> Flavors;
    typedef map<string, string> Plots;

    Flavors flavors;
    flavors.push_back("");
    flavors.push_back("b");
    flavors.push_back("c");
    flavors.push_back("uds");
    flavors.push_back("g");
    flavors.push_back("cl");
    flavors.push_back("l");

    Plots plots;
    plots["q"] ="other MuTag";
    plots["qtag"] ="other MuTag && Tagger";
    plots["n"] ="MuTag";
    plots["p"] ="MuTag && CMBtag";
    plots["ntag"] ="opp tag: MuTag";
    plots["ptag"] ="opp tag MuTag && CMBtag";
    plots["nnoTag"] ="opp tag: MuTag";
    plots["pnoTag"] ="opp tag MuTag && CMBtag";

    for(Flavors::const_iterator flavor = flavors.begin();
        flavors.end() != flavor;
        ++flavor)
    {
        fstore->add( new TH1D(("jet_deltaR" + (flavor->size() ? "_" + *flavor : "")).c_str(),"#Delta R",60,0.,0.55) );
        fstore->add( new TH1D(("jet_pTrel" + (flavor->size() ? "_" + *flavor : "")).c_str(),"p_{Trel} [GeV/c]" , 50, 0, 5 ) );
        fstore->add( new TH1F(("deltaPhi" + (flavor->size() ? "_" + *flavor : "")).c_str() ,"#Delta #phi",80,-3.15,3.15) );
	
        // Loop over different plot types: pT, Eta, etc.
        for(int i = 0; 3 > i; ++i)
        {
            string suffix;
            int nbins;
            Double_t *bins;

            switch(i)
            {
                case 0: suffix = "_pT";
                        nbins = nptbins;
                        bins = jetptbins;
                        break;

                case 1: suffix = "_eta";
                        nbins = netabins;
                        bins = jetetabins;
                        break;

                case 2: suffix = "_phi";
                        nbins = nphibins;
                        bins = jetphibins;
                        break;
            }

            for(Plots::const_iterator plot = plots.begin();
                plots.end() != plot;
                ++plot)
            {
                fstore->add( new TH2F((plot->first + suffix + (flavor->length() ? ("_" + *flavor) : "")).c_str(),
                                          (*flavor + suffix + plot->second).c_str(),
                                          nbins, bins, 50, 0., 5.),
                                      flavor->length() ? "MCTruth" : "muon_in_jet");
            }
        } // end loop over plot types
    } // End loop over flavors

    // Add All muon histograms
    //
    fstore->add( new TH1I( "all_muon_number", "All muons", 50, 0, 50) );
    fstore->add( new TH1F( "all_muon_pt", "All muon pt", 300, 0, 50) );
    fstore->add( new TH1F( "all_muon_phi", "All muon phi", 70, -3.5, 3.5) );
    fstore->add( new TH1F( "all_muon_eta", "All muon eta", 60, -3, 3) );

    // Add All jet histograms
    //
    fstore->add( new TH1I( "all_jet_number", "All jets", 50, 0, 50) );
    fstore->add( new TH1F( "all_jet_pt", "All jet pt", 300, 0, 50) );
    fstore->add( new TH1F( "all_jet_phi", "All jet phi", 70, -3.5, 3.5) );
    fstore->add( new TH1F( "all_jet_eta", "All jet eta", 60, -3, 3) );

/*
    fstore->add( new TH2F("n_pT","MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("p_pT","MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("ntag_pT","opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("ptag_pT","opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("nnoTag_pT","opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("pnoTag_pT","opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );

    fstore->add( new TH2F("q_pT","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("qtag_pT","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.), "muon_in_jet" );

    // Eta
    fstore->add( new TH2F("n_eta","MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("p_eta","MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );

    fstore->add( new TH2F("ntag_eta","opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("ptag_eta","opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("nnoTag_eta","opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("pnoTag_eta","opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );

    fstore->add( new TH2F("q_eta","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("qtag_eta","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.), "muon_in_jet" );

//
    fstore->add( new TH2F("n_pT_b","b MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_b","b MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_b","b opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_b","b opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_b","b opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_b","b opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_b","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_b","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_b","b MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_b","b MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_b","b opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_b","b opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_b","b opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_b","b opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_eta_b","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_b","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
//
    fstore->add( new TH2F("n_pT_cl","cl MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_cl","cl MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_cl","cl opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_cl","cl opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_cl","cl opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_cl","cl opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_cl","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_cl","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_cl","cl MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_cl","cl MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_cl","cl opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_cl","cl opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_cl","cl opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_cl","cl opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_eta_cl","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_cl","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
//
    fstore->add( new TH2F("n_pT_c","c MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_c","c MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_c","c opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_c","c opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_c","c opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_c","c opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_c","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_c","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_c","c MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_c","c MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_c","c opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_c","c opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_c","c opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_c","c opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_eta_c","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_c","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
//
    fstore->add( new TH2F("n_pT_l","l MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_l","l MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_l","l opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_l","l opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_l","l opp tag: MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_l","l opp tag MuTag && CMBtag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_l","other MuTag pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_l","other MuTag && Tagger pT vs pTrel",nptbins,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_l","l MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_l","l MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_l","l opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_l","l opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_l","l opp tag: MuTag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_l","l opp tag MuTag && CMBtag Eta vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("q_eta_l","other MuTag pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_l","other MuTag && Tagger pT vs pTrel",netabins,jetetabins,50,0.,5.), "MCTruth" );
    */

        /*
    fstore->add( new TH1D("jet_deltaR","#Delta R",60,0.,0.55) );
    fstore->add( new TH1D("jet_deltaR_b","#Delta R",60,0.,0.55) );
    fstore->add( new TH1D("jet_deltaR_c","#Delta R",60,0.,0.55) );
    */
    fstore->add( new TH1D("jet_deltaR_udsg","#Delta R",60,0.,0.55) );
	fstore->add( new TH1D("deltaRnearjet","#Delta R",60,0.,2.0) );
	
    // pending change colors
    /*
    fstore->add( new TH1D("jet_pTrel","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
    fstore->add( new TH1D("jet_pTrel_b","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
    fstore->add( new TH1D("jet_pTrel_c","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
    */
    fstore->add( new TH1D("jet_pTrel_udsg","p_{Trel} [GeV/c]" , 50, 0, 5 ) );

    fstore->add( new TH1F( "jet_pt", "jet pt", 30, 0, 150) );
    fstore->add( new TH1F( "awayjet_pt", "jet pt", 30, 0, 150) );
    fstore->add( new TH1F( "muon_pt", "muon pt", 300, 0, 50) );
    fstore->add( new TH1F( "ptRel", "ptRel", 100, 0, 10) );

    /*
    fstore->add( new TH1F( "deltaPhi","deltaPhi",80,-3.15,3.15) );
	fstore->add( new TH1F( "deltaPhi_b","deltaPhi_b",80,-3.15,3.15) );
	fstore->add( new TH1F( "deltaPhi_c","deltaPhi_c",80,-3.15,3.15) );
	fstore->add( new TH1F( "deltaPhi_l","deltaPhi_l",80,-3.15,3.15) );
	fstore->add( new TH1F( "deltaPhi_g","deltaPhi_g",80,-3.15,3.15) );
    */
}

//______________________________________________________________________________________________________________________
//void PMHistograms::FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
//			      int JetFlavor, std::map<std::string, bool> aMap)
void PMHistograms::FillHistos(const std::string &type, const TLorentzVector &p4MuJet, const double &ptrel,
                              int JetFlavor, bool tagged)
{
    // All jets
    //
    FillHisto("", type, p4MuJet, ptrel, tagged);

    // b jets
    //
    if ( JetFlavor == 5 )
        FillHisto("b", type, p4MuJet, ptrel, tagged);

    // c jets
    if ( JetFlavor == 4 )
        FillHisto("c", type, p4MuJet, ptrel, tagged);

    // uds jets
    //
    if ((JetFlavor > 0 && JetFlavor < 4))
        FillHisto("uds", type, p4MuJet, ptrel, tagged);

    // g jets
    if ( JetFlavor == 21 )
        FillHisto("g", type, p4MuJet, ptrel, tagged);

    // c + uds + g
    //
    if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
        FillHisto("cl", type, p4MuJet, ptrel, tagged);

    // uds + g
    //
    if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
        FillHisto("l", type, p4MuJet, ptrel, tagged);
}

void PMHistograms::FillHisto(const std::string &flavor,
                             const std::string &type,
                             const TLorentzVector &p4,
                             const double &ptrel,
                             const bool &tagged)
{
    string flvr = flavor.size()
                  ? "_" + flavor
                  : flavor;

    Fill(type, flvr, p4, ptrel);
    Fill(type + (tagged ? "tag" : "noTag"), flvr, p4, ptrel);
}

void PMHistograms::Fill(const std::string &prefix,
                        const std::string &suffix,
                        const TLorentzVector &p4,
                        const double &ptrel)
{
    fstore->hist(prefix + "_pT" + suffix)->Fill(p4.Pt(),ptrel);
    fstore->hist(prefix + "_eta" + suffix)->Fill(TMath::Abs(p4.Eta()), ptrel);
    fstore->hist(prefix + "_phi" + suffix)->Fill(TMath::Abs(p4.Phi()), ptrel);
}
