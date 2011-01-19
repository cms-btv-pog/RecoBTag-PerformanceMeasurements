/**
 * File, InputFile, OutputFile
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <memory>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

#include "IO/interface/FileDelegate.h"
#include "IO/interface/Event.h"
#include "S8Tree/interface/S8EventID.h"
#include "S8Tree/interface/S8Fwd.h"
#include "S8Tree/interface/S8GenEvent.h"

#include "IO/interface/File.h"

using std::clog;
using std::endl;

using s8::InputFile;
using s8::InputFileDelegate;
using s8::OutputFile;
using s8::OutputFileDelegate;

InputFile::InputFile() throw()
{
    _delegate = 0;
}

InputFile::~InputFile() throw()
{
    close();
}

InputFileDelegate *InputFile::delegate() const
{
    return _delegate;
}

void InputFile::setDelegate(InputFileDelegate *delegate)
{
    _delegate = delegate;
}

void InputFile::open()
{
    bool fileShouldOpen = true;
    if (_delegate)
        fileShouldOpen = _delegate->inputFileShouldOpen(name());

    if (!fileShouldOpen)
    {
        clog << "Skip file: " << name() << endl;

        return;
    }

    core::InputFile::open();

    if (_delegate)
        _delegate->inputFileDidOpen(file());
}

void InputFile::process()
{
    using std::cerr;

    clog << " Extract S8 Tree from input file" << endl;

    using std::runtime_error;

    TTree *tree = 0;
    file()->GetObject("S8TreeMaker/s8", tree);

    if (!tree)
        throw runtime_error("Failed to extract S8 Tree from: " + name());

    // All branches are enabled by default
    //

    std::auto_ptr<Event> event(new Event());

    // Create Branches
    //
    std::auto_ptr<EventID>         eventID(new EventID());
    std::auto_ptr<GenEvent>        genEvent(new GenEvent());
    std::auto_ptr<Jets>            jets(new Jets());
    std::auto_ptr<Leptons>         electrons(new Leptons());
    std::auto_ptr<Leptons>         muons(new Leptons());
    std::auto_ptr<PrimaryVertices> primaryVertices(new PrimaryVertices());
    std::auto_ptr<Triggers>        triggers(new Triggers());

    // Put branches into Event
    //
    event->setID(eventID.get());
    event->setGen(genEvent.get());
    event->setJets(jets.get());
    event->setElectrons(electrons.get());
    event->setMuons(muons.get());
    event->setPrimaryVertices(primaryVertices.get());
    event->setTriggers(triggers.get());

    // Enable branches
    //
    {
        s8::EventID *eventIDPointer = eventID.get();
        tree->SetBranchAddress("eventID", &eventIDPointer);
    }

    {
        s8::GenEvent *genEventPointer = genEvent.get();
        tree->SetBranchAddress("genEvent", &genEventPointer);
    }

    s8::Jets *jetsPointer = jets.get();
    tree->SetBranchAddress("jets", &jetsPointer);
    if (!_delegate ||
        !_delegate->inputFileShouldLoadJets())
    {
        clog << " + Disable Jets branch" << endl;

        tree->SetBranchStatus("*jets*", 0);
    }

    s8::Leptons *electronsPointer = electrons.get();
    tree->SetBranchAddress("electrons", &electronsPointer);
    if (!_delegate ||
        !_delegate->inputFileShouldLoadElectrons())
    {
        clog << " + Disable Electrons branch" << endl;

        tree->SetBranchStatus("*electrons*", 0);
    }

    s8::Leptons *muonsPointer = muons.get();
    tree->SetBranchAddress("muons", &muonsPointer);
    if (!_delegate ||
        !_delegate->inputFileShouldLoadMuons())
    {
        clog << " + Disable Muons branch" << endl;

        tree->SetBranchStatus("*muons*", 0);
    }

    s8::PrimaryVertices *primaryVerticesPointer = primaryVertices.get();
    tree->SetBranchAddress("primaryVertices", &primaryVerticesPointer);
    if (!_delegate ||
        !_delegate->inputFileShouldLoadPrimaryVertices())
    {
        clog << " + Disable PrimaryVertices branch" << endl;

        tree->SetBranchStatus("*primaryVertices*", 0);
    }

    s8::Triggers *triggersPointer = triggers.get();
    tree->SetBranchAddress("triggers", &triggersPointer);
    if (!_delegate ||
        !_delegate->inputFileShouldLoadTriggers())
    {
        clog << " + Disable Triggers branch" << endl;

        tree->SetBranchStatus("triggers", 0);
    }

    clog << endl;

    const int entries = tree->GetEntries();

    clog << " File has " << entries << " entries stored" << endl
        << endl;

    // Error counter
    //
    int errors = 0;

    // Loop over events
    //
    for(int entry = 0; entries > entry; ++entry)
    try
    {
        tree->GetEvent(entry);

        _delegate->inputFileDidLoadEvent(event.get());

        if (!_delegate->inputFileShouldContinue())
            break;
    }
    catch(...)
    {
        ++errors;

        if (5 > errors)
            continue;

        // Too many errors occured
        cerr << "5+ errors where fired while file is being read: exit" << endl;

        break;
    }

    clog << "done analyzing Tree" << endl;
}

void InputFile::close()
{
    if (_delegate)
        _delegate->inputFileWillClose(file());

    core::InputFile::close();
}



OutputFile::OutputFile() throw()
{
    _delegate = 0;
}

OutputFile::~OutputFile() throw()
{
    close();
}

OutputFileDelegate *OutputFile::delegate() const
{
    return _delegate;
}

void OutputFile::setDelegate(OutputFileDelegate *delegate)
{
    _delegate = delegate;
}

void OutputFile::open()
{
    core::OutputFile::open();

    if (_delegate)
        _delegate->outputFileDidOpen(file());
}

void OutputFile::close()
{
    if (_delegate)
        _delegate->outputFileWillClose(file());

    core::OutputFile::close();
}
