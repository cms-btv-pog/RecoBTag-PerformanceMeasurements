/**
 * MonitorAnalyzer + SolverInput
 * s8
 *
 * Created by Samvel Khalatian on Feb 14, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>

#include "Analyzer/interface/MonitorAnalyzer.h"
#include "Analyzer/interface/SolverInputAnalyzer.h"

#include "Analyzer/interface/MonitorSolverInputAnalyzer.h"

using std::cout;
using std::endl;

using s8::MonitorSolverInputAnalyzer;

MonitorSolverInputAnalyzer::MonitorSolverInputAnalyzer() throw()
{
    _monitor.reset(new MonitorAnalyzer());
    _solver.reset(new SolverInputAnalyzer());
}

MonitorSolverInputAnalyzer::~MonitorSolverInputAnalyzer() throw()
{
}



// Analyzer interface
//
void MonitorSolverInputAnalyzer::init()
{
    _monitor->init();
    _solver->init();
}

void MonitorSolverInputAnalyzer::treeDidLoad(const TreeInfo *tree_info,
                                             const TriggerCenter *trigger_center)
{
    _monitor->treeDidLoad(tree_info, trigger_center);
    _solver->treeDidLoad(tree_info, trigger_center);
}

void MonitorSolverInputAnalyzer::eventDidLoad(const Event *event)
{
    _monitor->eventDidLoad(event);
    _solver->eventDidLoad(event);
}

void MonitorSolverInputAnalyzer::print(std::ostream &out) const
{
    _monitor->print(out);
    _solver->print(out);
}

void MonitorSolverInputAnalyzer::save(TDirectory *directory) const
{
    _monitor->save(directory);
    _solver->save(directory);
}

// SolverInputOptionsDelegate interface
//
void MonitorSolverInputAnalyzer::optionDataIsSet(const bool &value)
{
    _solver->optionDataIsSet(value);
}

// MuonInJetOptionsDelegate interface
//
void MonitorSolverInputAnalyzer::optionTagIsSet(const std::string &tag)
{
    _monitor->optionTagIsSet(tag);
    _solver->optionTagIsSet(tag);
}

void MonitorSolverInputAnalyzer::optionAwayTagIsSet(const std::string &tag)
{
    _monitor->optionAwayTagIsSet(tag);
    _solver->optionAwayTagIsSet(tag);
}

void MonitorSolverInputAnalyzer::optionMuonPtIsSet(const Range &value)
{
    _monitor->optionMuonPtIsSet(value);
    _solver->optionMuonPtIsSet(value);
}

void MonitorSolverInputAnalyzer::optionJetPtIsSet(const Range &value)
{
    _monitor->optionJetPtIsSet(value);
    _solver->optionJetPtIsSet(value);
}

void MonitorSolverInputAnalyzer::optionJetEtaIsSet(const Range &value)
{
    _monitor->optionJetEtaIsSet(value);
    _solver->optionJetEtaIsSet(value);
}

// PythiaOptionsDelegate interface
//
void MonitorSolverInputAnalyzer::optionGluonSplittingIsSet(const GluonSplitting &value)
{
    _monitor->optionGluonSplittingIsSet(value);
    _solver->optionGluonSplittingIsSet(value);
}

void MonitorSolverInputAnalyzer::optionPtHatIsSet(const Range &value)
{
    _monitor->optionPtHatIsSet(value);
    _solver->optionPtHatIsSet(value);
}

// Trigger options
//
void MonitorSolverInputAnalyzer::optionTriggerIsSet(const Trigger &trigger)
{
    _monitor->optionTriggerIsSet(trigger);
    _solver->optionTriggerIsSet(trigger);
}

void MonitorSolverInputAnalyzer::optionSimulateTriggerIsSet(const bool &value)
{
    _monitor->optionSimulateTriggerIsSet(value);
    _solver->optionSimulateTriggerIsSet(value);
}

void MonitorSolverInputAnalyzer::optionReweightTriggerIsSet(const std::string &filename)
{
    _monitor->optionReweightTriggerIsSet(filename);
    _solver->optionReweightTriggerIsSet(filename);
}

// Misc Options
//
void MonitorSolverInputAnalyzer::optionPrimaryVerticesIsSet(const Range &primary_vertices)
{
    _monitor->optionPrimaryVerticesIsSet(primary_vertices);
    _solver->optionPrimaryVerticesIsSet(primary_vertices);
}
