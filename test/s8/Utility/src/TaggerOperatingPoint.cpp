/**
 * TaggerOperatingPoint
 * s8
 *
 * Created by Samvel Khalatian on Nov 15, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include "Utility/interface/TaggerOperatingPoint.h"

using std::runtime_error;

using s8::TaggerOperatingPoint;

TaggerOperatingPoint::TaggerOperatingPoint() throw()
{
    initWithOperatingPoint(TCHPL);
}

TaggerOperatingPoint::~TaggerOperatingPoint() throw()
{
}

void TaggerOperatingPoint::initWithOperatingPoint(const OperatingPoint &op)
{
    switch(op)
    {
        case TCHET: _operatingPoint = 10.2;
                    _btag = Jet::TCHE;
                    break;

        case TCHEM: _operatingPoint = 3.3;
                    _btag = Jet::TCHE;
                    break;

        case TCHEL: _operatingPoint = 1.7;
                    _btag = Jet::TCHE;
                    break;



        case TCHPT: _operatingPoint = 3.41;
                    _btag = Jet::TCHP;
                    break;

        case TCHPM: _operatingPoint = 1.93;
                    _btag = Jet::TCHP;
                    break;

        case TCHPL: _operatingPoint = 1.19;
                    _btag = Jet::TCHP;
                    break;



        case SSVT:  _operatingPoint = 3.05;
                    _btag = Jet::SSV;
                    break;

        case SSVM:  _operatingPoint = 1.74;
                    _btag = Jet::SSV;
                    break;

        case SSVL:  _operatingPoint = 1.05;
                    _btag = Jet::SSV;
                    break;



        case SSVHET: _operatingPoint = 3.05;
                     _btag = Jet::SSVHE;
                     break;

        case SSVHEM: _operatingPoint = 1.74;
                     _btag = Jet::SSVHE;
                     break;



        case SSVHPT: _operatingPoint = 2.00;
                     _btag = Jet::SSVHP;
                     break;



        default:
                     throw runtime_error("Unsupported OperatingPoint supplied");
    }
}

TaggerOperatingPoint::operator double() const
{
    return _operatingPoint;
}

s8::Jet::BTag TaggerOperatingPoint::btag() const
{
    return _btag;
}



// Helpers
//
TaggerOperatingPoint &s8::operator<<(TaggerOperatingPoint &op,
                                     const std::string &value)
{
    if ("TCHET" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::TCHET);

    else if ("TCHEM" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::TCHEM);

    else if ("TCHEL" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::TCHEL);



    else if ("TCHPT" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::TCHPT);

    else if ("TCHPM" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::TCHPM);

    else if ("TCHPL" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::TCHPL);



    else if ("SSVT" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::SSVT);

    else if ("SSVM" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::SSVM);

    else if ("SSVL" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::SSVL);



    else if ("SSVHET" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::SSVHET);

    else if ("SSVHEM" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::SSVHEM);



    else if ("SSVHPT" == value)
        op.initWithOperatingPoint(TaggerOperatingPoint::SSVHPT);


    else
        throw runtime_error("Unknown Tagger Operating Point defined: " + value);

    return op;
}

bool s8::operator==(const TaggerOperatingPoint::OperatingPoint &op1,
                    const TaggerOperatingPoint &op2)
{
    TaggerOperatingPoint op;
    op.initWithOperatingPoint(op1);

    return op == op2;
}

bool s8::operator==(const TaggerOperatingPoint &op1,
                    const TaggerOperatingPoint &op2)
{
    return op1.btag() == op2.btag() &&
           op1.operator double() == op2.operator double();
}

bool s8::operator!=(const TaggerOperatingPoint &op1,
                    const TaggerOperatingPoint &op2)
{
    return !(op1 == op2);
}
