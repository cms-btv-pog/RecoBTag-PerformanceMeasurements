#include "RecoBTag/PerformanceMeasurements/interface/AnalyzerBase.h"

bool analyzerBase::compare(const P& i, const P& j) {
    return i.second.index() > j.second.index();
}

const reco::TrackBaseRef analyzerBase::toTrackRef(const reco::TrackRef & trk) {return reco::TrackBaseRef(trk);}
const reco::TrackBaseRef analyzerBase::toTrackRef(const edm::Ptr<reco::Candidate> & cnd)
{
  const pat::PackedCandidate * pcand = dynamic_cast<const pat::PackedCandidate *>(cnd.get());

  if(pcand) // MiniAOD case
    return reco::TrackBaseRef(); // return null reference since no tracks are stored in MiniAOD
  else
  {
    const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(cnd.get());

    if ( (std::abs(pfcand->pdgId()) == 11 || pfcand->pdgId() == 22) && pfcand->gsfTrackRef().isNonnull() && pfcand->gsfTrackRef().isAvailable() )
      return reco::TrackBaseRef(pfcand->gsfTrackRef());
    else if ( pfcand->trackRef().isNonnull() && pfcand->trackRef().isAvailable() )
      return reco::TrackBaseRef(pfcand->trackRef());
    else
      return reco::TrackBaseRef();
  }
}

const math::XYZPoint & analyzerBase::position(const reco::Vertex & sv) {return sv.position();}
const math::XYZPoint & analyzerBase::position(const reco::VertexCompositePtrCandidate & sv) {return sv.vertex();}
const double analyzerBase::xError(const reco::Vertex & sv) {return sv.xError();}
const double analyzerBase::xError(const reco::VertexCompositePtrCandidate & sv) {return sqrt(sv.vertexCovariance(0,0));}
const double analyzerBase::yError(const reco::Vertex & sv) {return sv.yError();}
const double analyzerBase::yError(const reco::VertexCompositePtrCandidate & sv) {return sqrt(sv.vertexCovariance(1,1));}
const double analyzerBase::zError(const reco::Vertex & sv) {return sv.zError();}
const double analyzerBase::zError(const reco::VertexCompositePtrCandidate & sv) {return sqrt(sv.vertexCovariance(2,2));}
const double analyzerBase::chi2(const reco::Vertex & sv) {return sv.chi2();}
const double analyzerBase::chi2(const reco::VertexCompositePtrCandidate & sv) {return sv.vertexChi2();}
const double analyzerBase::ndof(const reco::Vertex & sv) {return sv.ndof();}
const double analyzerBase::ndof(const reco::VertexCompositePtrCandidate & sv) {return sv.vertexNdof();}
const unsigned int analyzerBase::vtxTracks(const reco::Vertex & sv) {return sv.nTracks();}
const unsigned int analyzerBase::vtxTracks(const reco::VertexCompositePtrCandidate & sv) {return sv.numberOfSourceCandidatePtrs();}


analyzerBase::JetFlavor analyzerBase::jet_flavour(const pat::Jet& jet,
    const std::vector<reco::GenParticle>& gToBB,
    const std::vector<reco::GenParticle>& gToCC,
    const std::vector<reco::GenParticle>& neutrinosLepB,
    const std::vector<reco::GenParticle>& neutrinosLepB_C,
    const std::vector<reco::GenParticle>& alltaus,
    bool usePhysForLightAndUndefined) {
    int hflav = abs(jet.hadronFlavour());
    int pflav = abs(jet.partonFlavour());
    int physflav = 0;
    if(jet.genParton()) physflav=abs(jet.genParton()->pdgId());
    std::size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
    std::size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

    unsigned int nbFromGSP(0);
    for (reco::GenParticle p : gToBB) {
        double dr(reco::deltaR(jet, p));
        if (dr < 0.4) ++nbFromGSP;
    }

    unsigned int ncFromGSP(0);
    for (reco::GenParticle p : gToCC) {
        double dr(reco::deltaR(jet, p));
        if (dr < 0.4) ++ncFromGSP;
    }

    //std::cout << " jet pt = " << jet.pt() << " hfl = " << hflav << " pfl = " << pflav << " genpart = " << physflav
            //  << " nbFromGSP = " << nbFromGSP << " ncFromGSP = " << ncFromGSP
    //  << " nBhadrons " << nbs << " nCHadrons " << ncs << std::endl;

    if(hflav == 5) { //B jet
        if(nbs > 1) {
            if (nbFromGSP > 0) return analyzerBase::JetFlavor::GBB;
            else return analyzerBase::JetFlavor::BB;
        }
        else if(nbs == 1) {
            for (std::vector<reco::GenParticle>::const_iterator it = neutrinosLepB.begin(); it != neutrinosLepB.end(); ++it){
                if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4) {
                    return analyzerBase::JetFlavor::LeptonicB;
                }
            }
            for (std::vector<reco::GenParticle>::const_iterator it = neutrinosLepB_C.begin(); it != neutrinosLepB_C.end(); ++it){
                if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4) {
                    return analyzerBase::JetFlavor::LeptonicB_C;
                }
            }
            return analyzerBase::JetFlavor::B;
        }
        else {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return analyzerBase::JetFlavor::G;
                else if(physflav == 3) return analyzerBase::JetFlavor::S;
                else if(physflav == 2 || physflav ==1) return analyzerBase::JetFlavor::UD;
                else return analyzerBase::JetFlavor::UNDEFINED;
            }
            else return analyzerBase::JetFlavor::UNDEFINED;
        }
    }
    else if(hflav == 4) { //C jet
        if (ncs > 1) {
            if (ncFromGSP > 0) return analyzerBase::JetFlavor::GCC;
            else return analyzerBase::JetFlavor::CC;
        }
        else return analyzerBase::JetFlavor::C;
    }
    else { //not a heavy jet
        if(alltaus.size()>0){ //check for tau in a simplistic way
            bool ishadrtaucontained=true;
            for(const auto& p:alltaus){
                size_t ndau=p.numberOfDaughters();
                for(size_t i=0;i<ndau;i++){
                    const reco::Candidate* dau=p.daughter(i);
                    int daupid=std::abs(dau->pdgId());
                    if(daupid == 13 || daupid == 11){
                        ishadrtaucontained=false;
                        break;
                    }
                    if(daupid != 12 && daupid!=14 && daupid!=16 &&
                            reco::deltaR(*dau,jet) > 0.4){
                        ishadrtaucontained=false;
                        break;
                    }
                }
            }
            if(ishadrtaucontained) return analyzerBase::JetFlavor::TAU;
        }
        if(std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs) {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return analyzerBase::JetFlavor::G;
                else if(physflav == 3) return analyzerBase::JetFlavor::S;
                else if(physflav == 2 || physflav ==1) return analyzerBase::JetFlavor::UD;
                else return analyzerBase::JetFlavor::UNDEFINED;
            }
            else return analyzerBase::JetFlavor::UNDEFINED;
        }
        else if(usePhysForLightAndUndefined){
            if(physflav == 21) return analyzerBase::JetFlavor::G;
            else if(physflav == 3) return analyzerBase::JetFlavor::S;
            else if(physflav == 2 || physflav ==1) return analyzerBase::JetFlavor::UD;
            else return analyzerBase::JetFlavor::UNDEFINED;
        }
        else {
            if(pflav == 21) return analyzerBase::JetFlavor::G;
            else if(pflav == 3) return analyzerBase::JetFlavor::S;
            else if(pflav == 2 || pflav ==1) return analyzerBase::JetFlavor::UD;
            else return analyzerBase::JetFlavor::UNDEFINED;
        }
    }
    return analyzerBase::JetFlavor::UNDEFINED;
}
