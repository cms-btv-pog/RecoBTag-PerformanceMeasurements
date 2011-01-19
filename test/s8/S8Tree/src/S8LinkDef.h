#include "interface/S8GenEvent.h"
#include "interface/S8GenParticle.h"
#include "interface/S8EventID.h"
#include "interface/S8Jet.h"
#include "interface/S8Lepton.h"
#include "interface/S8PrimaryVertex.h"
#include "interface/S8TreeInfo.h"
#include "interface/S8Trigger.h"
#include "interface/S8TriggerCenter.h"
#include "interface/S8Fwd.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace s8;

#pragma link C++ class s8::EventID+;
#pragma link C++ class s8::GenEvent+;
#pragma link C++ class s8::GenParticle+;
#pragma link C++ class s8::Lepton+;
#pragma link C++ class s8::Jet+;
#pragma link C++ class s8::PrimaryVertex+;
#pragma link C++ class s8::Trigger+;

#pragma link C++ class std::vector<s8::Jet *>;
#pragma link C++ class std::vector<s8::Lepton *>;
#pragma link C++ class std::vector<s8::PrimaryVertex *>;
#pragma link C++ class std::vector<s8::Trigger *>;

#pragma link C++ class s8::Version+;
#pragma link C++ class s8::TreeInfo+;

#pragma link C++ class std::pair<unsigned int, std::string>;
#pragma link C++ class std::map<unsigned int, std::string>;
#pragma link C++ class s8::TriggerCenter+;

#endif
