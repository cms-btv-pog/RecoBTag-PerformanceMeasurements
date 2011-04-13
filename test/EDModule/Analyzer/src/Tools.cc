/**
 * Tools
 * top
 *
 * Created by Samvel Khalatian on Sep 7, 2010
 * Copyright 2010, All rights reserved
 */

#include <TLorentzVector.h>
#include <TVector3.h>

#include "DataFormats/MuonReco/interface/Muon.h"

#include "Tree/Top/interface/TopJet.h"
#include "Tree/Top/interface/TopJetEnergy.h"
#include "Tree/Top/interface/TopElectronIsolation.h"
#include "Tree/Top/interface/TopMuonIsolation.h"

#include "EDModule/Analyzer/interface/Tools.h"

namespace tt = top::tools;

void tt::setP4(TLorentzVector *topP4,
               const math::XYZTLorentzVector *cmsswP4)
{
    topP4->SetPxPyPzE(cmsswP4->px(),
                      cmsswP4->py(),
                      cmsswP4->pz(),
                      cmsswP4->energy());
}

void tt::setVertex(TVector3 &vector, const math::XYZPoint &point)
{
    vector.SetXYZ(point.x(), point.y(), point.z());
}

void tt::setEnergy(top::Jet &jet,
                   const reco::CaloJet::Specific &specific)
{
    top::JetEnergy energy;

    energy.setEcalMax(specific.mMaxEInEmTowers);
    energy.setHcalMax(specific.mMaxEInHadTowers);

    energy.setHcalInHO(specific.mHadEnergyInHO);
    energy.setHcalInHB(specific.mHadEnergyInHB);
    energy.setHcalInHF(specific.mHadEnergyInHF);
    energy.setHcalInHE(specific.mHadEnergyInHE);

    energy.setEcalInEB(specific.mEmEnergyInEB);
    energy.setEcalInEE(specific.mEmEnergyInEE);
    energy.setEcalInHF(specific.mEmEnergyInHF);

    energy.setEcalFraction(specific.mEnergyFractionEm);
    energy.setHcalFraction(specific.mEnergyFractionHadronic);

    jet.setEnergy(energy);
}

void tt::setIsolation(top::Muon &muon,
                      const top::Muon::ISO &iso,
                      const reco::MuonIsolation &recoIso)
{
    top::MuonIsolation topIso;

    topIso.setTrackPt(recoIso.sumPt);
    topIso.setEcalEt(recoIso.emEt);
    topIso.setHcalEt(recoIso.hadEt);

    topIso.setTracks(recoIso.nTracks);
    topIso.setJets(recoIso.nJets);

    topIso.setTrackPtVeto(recoIso.trackerVetoPt);
    topIso.setEcalEtVeto(recoIso.emVetoEt);
    topIso.setHcalEtVeto(recoIso.hadVetoEt);

    muon.setIsolation(iso, topIso);
}

void tt::setIsolation(top::Electron &electron,
                      const top::Electron::ISO &iso,
                      const reco::GsfElectron::IsolationVariables &recoIso)
{
    top::ElectronIsolation topIso;

    topIso.setTrackPt(recoIso.tkSumPt);
    topIso.setEcalEt(recoIso.ecalRecHitSumEt);
    topIso.setHcalEt(recoIso.hcalDepth1TowerSumEt + recoIso.hcalDepth2TowerSumEt);

    electron.setIsolation(iso, topIso);
}
