/**
 * File, InputFile, OutputFile
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include <TFile.h>

#include "interface/File.h"

using std::runtime_error;
using std::string;

using core::File;
using core::InputFile;
using core::OutputFile;

File::File() throw()
{
}

File::~File() throw()
{
}

string File::name() const
{
    return _name;
}

void File::setName(const string &name)
{
    _name = name;
}

void File::open(const Mode &mode)
{
    if (_name.empty())
        throw runtime_error("Input filename is not set");

    if (mode & RECREATE)
        _file.reset(new TFile(_name.c_str(), "recreate"));
    else
        _file.reset(new TFile(_name.c_str()));

    if (!_file.get() ||
        !_file->IsOpen())
    {
        _file.reset();

        throw runtime_error("Failed to open input file: " + _name);
    }
}

bool File::isOpen() const
{
    return file() ? true : false;
}

void File::close()
{
    _file.reset();
}

TFile *File::file() const
{
    return _file.get();
}



InputFile::InputFile() throw()
{
}

InputFile::~InputFile() throw()
{
}

void InputFile::open()
{
    File::open(File::READ);
}



OutputFile::OutputFile() throw()
{
}

OutputFile::~OutputFile() throw()
{
}

void OutputFile::open()
{
    File::open(File::RECREATE);
}
