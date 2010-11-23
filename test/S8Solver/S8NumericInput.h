/**
 * S8NumericInput
 * 
 *
 * Created by Samvel Khalatian on Nov 10, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_NUMERICINPUT
#define S8_NUMERICINPUT

#include <ostream>
#include <vector>
#include <utility>

// Central value and variance (binomial)
//
typedef std::pair<double, double> Measurement;

namespace numeric
{
    struct InputGroup
    {
        Measurement all;
        Measurement tag;
        Measurement mu;
        Measurement muTag;
    };



    struct FlavouredInput
    {
        virtual ~FlavouredInput() {}

        Measurement b;
        Measurement cl;
    };

    struct FlavouredInputGroup: public FlavouredInput
    {
        FlavouredInput tag;
        FlavouredInput mu;
        FlavouredInput muTag;
    };
}

struct NumericInput
{
    numeric::InputGroup n;
    numeric::InputGroup p;
};

struct FlavouredNumericInput
{
    numeric::FlavouredInputGroup n;
    numeric::FlavouredInputGroup p;
};

struct FlavouredEffGroup
{
    numeric::FlavouredInput tag;
    numeric::FlavouredInput mu;
};

struct Coefficients
{
    Measurement alpha;
    Measurement beta;
    Measurement gamma;
    Measurement delta;
    Measurement kappaCL;
    Measurement kappaB;
};

struct NumericInputGroup
{
    NumericInput          input;
    FlavouredNumericInput flavouredInput;
    FlavouredEffGroup     efficiency;
    Coefficients          coefficients;
    Measurement           bin;
};

std::ostream &operator<<(std::ostream &, const Measurement &);

Measurement &operator+=(Measurement &, const Measurement &);

Measurement operator +(const Measurement &, const Measurement &);

Measurement operator *(const double &, const Measurement &);

// Stat error propagation: indep measurements
//
Measurement operator /(const Measurement &,
                       const Measurement &);

// Binomial error calculation
//
Measurement operator %(const Measurement &,
                       const Measurement &);


std::ostream &operator<<(std::ostream &, const NumericInput &);

std::ostream &operator<<(std::ostream &, const FlavouredNumericInput &);

std::ostream &operator<<(std::ostream &, const FlavouredEffGroup &);

std::ostream &operator<<(std::ostream &, const Coefficients &);

std::ostream &operator<<(std::ostream &, const NumericInputGroup &);

#endif
