
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "Math/GenVector/VectorUtil.h"

#include "RecoBTag/PerformanceMeasurements/interface/PFTools.h"


namespace PFTools
{


reco::GenJet GetGenJet(reco::CaloJet calojet, reco::GenJetCollection genJetColl)
{
    reco::GenJet matchedJet;

    double predelta = 99999.;
    for (reco::GenJetCollection::const_iterator genjet = genJetColl.begin(); genjet != genJetColl.end(); ++genjet)
    {
        double delta  = ROOT::Math::VectorUtil::DeltaR(genjet->p4().Vect(), calojet.p4().Vect() );

        if ( delta < 0.2 && delta<predelta )
        {
            matchedJet = *genjet;
            predelta = delta;
        }
    }

    return matchedJet;
}


int TaggedJet(reco::CaloJet const & calojet, edm::Handle<reco::JetTagCollection> const & jetTags)
{

    double small = 1.e-5;
    int result = -1; // no tagged
    //int ith = -1;

    //std::cout << "calo jet: pz = " << calojet.pz() << " pt = " << calojet.pt() << std::endl;
    //for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
    //  edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];

    //get label and module names


    //    std::cout <<" ECCO " << jetTags.product()<< std::endl;


    for (size_t t = 0; t < jetTags->size(); ++t)
    {
        edm::RefToBase<reco::Jet> jet_p = (*jetTags)[t].first;
        if (jet_p.isNull())
        {
            //std::cout << "-----------> JetTag::jet() returned null reference" << std::endl;
            continue;
        }
        //std::cout << "[TaggedJet]  calojet pt = " << calojet.pt() << " tagged jet pt = " << jet_p->pt() << std::endl;
        if (DeltaR<reco::Candidate>()( calojet, *jet_p ) < small)
        {

            result = (int) t;

        }
    }

    return result;
}


int TaggedJet(reco::CaloJet const & caloJet, edm::Handle<std::vector<reco::TrackIPTagInfo> > const & trackIPTagInfos )
{
    double small = 1.e-5;
    int result = -1;

    for (size_t t = 0; t < trackIPTagInfos->size(); ++t)
    {
        edm::RefToBase<reco::Jet> jet_p = (*trackIPTagInfos)[t].jet();
        if (jet_p.isNull())
        {
            //std::cout << "-----------> TrackIPTagInfos::jet() returned null reference" << std::endl;
            continue;
        }
        //std::cout << "[TaggedJet]  calojet pt = " << calojet.pt() << " tagged jet pt = " << jet_p->pt() << std::endl;
        if (DeltaR<reco::Candidate>()( caloJet, * jet_p ) < small)
        {
            result = (int) t;
            break;
        }
    }

    return result;
}


std::map<std::string, bool> GetBTaggingMap(std::vector<WorkingPoint> const & wp, reco::CaloJet const & jet, edm::Event const & event, double ptrel)
{
    std::map<std::string, bool> aMap;

    int ith_tagged = -1;

    for (std::vector<WorkingPoint>::const_iterator it = wp.begin(); it != wp.end(); ++it)
    {

        edm::Handle<reco::JetTagCollection > jetTags;
        event.getByLabel((*it).inputTag(),jetTags);

        //std::string moduleLabel = (jetTags).provenance()->moduleLabel();
        //if (mymap.find(moduleLabel) != mymap.end()) continue;
        //mymap[moduleLabel] = true;

        ith_tagged = TaggedJet(jet,jetTags);

        if (ith_tagged == -1) continue;

        std::map<std::string, double > list_cuts = (*it).list();

        for (std::map<std::string, double >::const_iterator icut = list_cuts.begin(); icut != list_cuts.end(); ++icut)
        {
            std::string alabel = icut->first;
            if ( (*jetTags)[ith_tagged].second > icut->second ) aMap[alabel] = true;
            else aMap[alabel] = false;

        }

        //here
        /*
        edm::Handle<reco::JetTagCollection> jetTags;
        event.getByLabel((*it).inputTag(), jetTags);


        ith_tagged = TaggedJet(jet, jetTags);

        if (ith_tagged == -1) continue;

        if ((*jetTags)[ith_tagged].second > (*it).cut())
        {
            aMap[(*it).name()] = true;
        }
        else
        {
            aMap[(*it).name()] = false;
        }
        */

    }

    return aMap;
}



// THis function seems to not be use, if that the case please remove it !!!
/* SimTrack GetGenTrk(reco::Track const & atrack, const edm::SimTrackContainer *simTrkColl, const edm::SimVertexContainer *simVtcs)
{
    SimTrack matchedTrk;
    edm::SimVertexContainer mysimVtcs = *simVtcs;

    double predelta = 99999.;
    for (edm::SimTrackContainer::const_iterator gentrk = simTrkColl->begin(); gentrk != simTrkColl->end(); ++gentrk)
    {
        double delta = reco::deltaR( (*gentrk).momentum(), atrack.momentum() );

        int type = (*gentrk).type();
        if ( delta < 0.2 && delta<predelta && ((*gentrk).charge() == atrack.charge() ) &&
                ( abs(type)==11 || abs(type)==13 || abs(type)==15 || abs(type)==211 || abs(type)==321 ) )
        {
            matchedTrk = *gentrk;
            predelta = delta;
            math::XYZTLorentzVectorD v = (mysimVtcs)[(*gentrk).vertIndex()].position();

            //std::cout << "gentrk: vx = " << v.x() << std::endl;
            //std::cout << "rectrk: vx = " << atrack.vx() << std::endl;

        }
    }

    return matchedTrk;
}

int PerformanceAnalyzer::GetMotherId(const edm::SimVertexContainer *simVtxColl, const edm::SimTrackContainer *simTrkColl, SimTrack muonMC)
{

    edm::SimVertexContainer mysimVtxColl = *simVtxColl;

    edm::SimTrackContainer mysimTrkColl = *simTrkColl;

    // fill map of simtracks
    std::map<unsigned, unsigned> geantToIndex;
    for ( unsigned it=0; it< mysimTrkColl.size(); ++it )
    {
        geantToIndex[ mysimTrkColl[it].trackId() ] = it;
    }

    // The origin vertex
    int vertexId = muonMC.vertIndex();
    SimVertex vertex = mysimVtxColl[vertexId];

    // The mother track
    int motherId = -1;
    if ( vertex.parentIndex() )  // there is a parent to this vertex
    {

        // geant id of the mother
        unsigned motherGeandId =   vertex.parentIndex();
        std::map<unsigned, unsigned >::iterator association
        = geantToIndex.find( motherGeandId );
        if (association != geantToIndex.end() )
            motherId = association->second;
    }


    int motherType = motherId == -1 ? 0 : mysimTrkColl[motherId].type();

    return motherType;

}*/

}
