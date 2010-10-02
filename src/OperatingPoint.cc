/**
 * Operating Point
 * 
 *
 * Created by Samvel Khalatian on Oct 2, 2010
 * Copyright 2010, All rights reserved
 */

#include <stdexcept>

#include "RecoBTag/PerformanceMeasurements/interface/OperatingPoint.h"

using std::runtime_error;
using std::string;

using s8::OperatingPoint;

OperatingPoint::OperatingPoint()
{
    set(TCHEM);
}

void OperatingPoint::set(const OperatingPoint::OP &op)
{
    switch(op)
    {
        case TCHET:  _value = 10.2;
                     break;

        case TCHEM:  _value = 3.3;
                     break;

        case TCHEL:  _value = 1.7;
                     break;



        case TCHPT:  _value = 3.41;
                     break;

        case TCHPM:  _value = 1.93;
                     break;

        case TCHPL:  _value = 1.19;
                     break;



        case SSVT:   _value = 3.05;
                     break;

        case SSVM:   _value = 1.74;
                     break;

        case SSVL:   _value = 1.05;
                     break;



        case SSVHET: _value = 3.05;
                     break;

        case SSVHEM: _value = 1.74;
                     break;



        case SSVHPT: _value = 2.00;
                     break;



        default:     throw runtime_error("Bad OperatingPoint supplied.");
    }
}

OperatingPoint &s8::operator<<(OperatingPoint &op, const string &opString)
{
    if ("TCHET" == opString)
        op.set(OperatingPoint::TCHET);

    else if ("TCHEM" == opString)
        op.set(OperatingPoint::TCHEM);

    else if ("TCHEL" == opString)
        op.set(OperatingPoint::TCHEL);



    else if ("TCHPT" == opString)
        op.set(OperatingPoint::TCHPT);

    else if ("TCHPM" == opString)
        op.set(OperatingPoint::TCHPM);

    else if ("TCHPL" == opString)
        op.set(OperatingPoint::TCHPL);



    else if ("SSVT" == opString)
        op.set(OperatingPoint::SSVT);

    else if ("SSVM" == opString)
        op.set(OperatingPoint::SSVM);

    else if ("SSVL" == opString)
        op.set(OperatingPoint::SSVL);



    else if ("SSVHET" == opString)
        op.set(OperatingPoint::SSVHET);

    else if ("SSVHEM" == opString)
        op.set(OperatingPoint::SSVHEM);



    else if ("SSVHPT" == opString)
        op.set(OperatingPoint::SSVHPT);


    else
        throw runtime_error("Unknown operating point defined: " + opString);

    return op;
}
