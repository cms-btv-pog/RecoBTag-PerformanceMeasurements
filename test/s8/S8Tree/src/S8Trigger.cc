/**
 * Trigger
 * s8 
 *
 * Created by Samvel Khalatian on Oct 22, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <stdexcept>

#include "interface/S8Trigger.h"

using std::runtime_error;

using s8::Trigger;

Trigger::Trigger() throw():
    _hash(0),
    _version(1),
    _prescale(1),
    _isPass(false)
{
}

Trigger::Hash Trigger::hash() const
{
    return _hash;
}

int Trigger::version() const
{
    return _version;
}

int Trigger::prescale() const
{
    return _prescale;
}

Trigger::operator bool() const
{
    return _isPass;
}

void Trigger::setHash(const Hash &hash)
{
    _hash = hash;
}

void Trigger::setVersion(const int &version)
{
    if (0 == _hash)
        throw runtime_error("Trigger is undefined");

    if (1 > version)
        throw runtime_error("Version can not be negative");

    _version = version;
}

void Trigger::setPrescale(const int &prescale)
{
    if (0 > prescale)
        throw runtime_error("Negative Prescale supplied");

    _prescale = prescale;
}

void Trigger::setIsPass(const bool &isPass)
{
    if (0 == _hash)
        throw runtime_error("Trigger is undefined");

    _isPass = isPass;
}
