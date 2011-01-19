/**
 * InputFileDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include "IO/interface/FileDelegate.h"

using s8::InputFileDelegate;
using s8::OutputFileDelegate;

InputFileDelegate::~InputFileDelegate() throw()
{
}

// Always open file by default
//
bool InputFileDelegate::inputFileShouldOpen(const std::string &)
{
    return true;
}

void InputFileDelegate::inputFileDidOpen(TFile *)
{
}

bool InputFileDelegate::inputFileShouldLoadJets()
{
	return false;
}

bool InputFileDelegate::inputFileShouldLoadElectrons()
{
	return false;
}

bool InputFileDelegate::inputFileShouldLoadMuons()
{
	return false;
}

bool InputFileDelegate::inputFileShouldLoadPrimaryVertices()
{
	return false;
}

bool InputFileDelegate::inputFileShouldLoadTriggers()
{
	return false;
}

void InputFileDelegate::inputFileDidLoadEvent(const Event *)
{
}

bool InputFileDelegate::inputFileShouldContinue()
{
    return true;
}

void InputFileDelegate::inputFileWillClose(TFile *)
{
}



OutputFileDelegate::~OutputFileDelegate() throw()
{
}

void OutputFileDelegate::outputFileDidOpen(TFile *)
{
}

void OutputFileDelegate::outputFileWillClose(TFile *)
{
}
