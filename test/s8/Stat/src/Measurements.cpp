/**
 * Measurement
 * s8
 *
 * Created by Samvel Khalatian on Dec 20, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>

#include <ostream>
#include <stdexcept>

#include "Stat/interface/Measurement.h"

using std::runtime_error;

using s8::Measurement;

Measurement::Measurement() throw():
    _value(0),
    _variance(0)
{
}

Measurement::Measurement(const double &value, const double &variance) throw():
    _value(value),
    _variance(variance)
{
}

double Measurement::value() const
{
    return _value;
}

double Measurement::variance() const
{
    return _variance;
}

void Measurement::setValue(const double &value)
{
    _value = value;
}

void Measurement::setVariance(const double &variance)
{
    _variance = variance;
}



// Helpers
//
std::ostream &s8::operator<<(std::ostream &out, const Measurement &m)
{
    using std::ios_base;
    using std::ios;

    const std::streamsize precision = out.precision(4);
    const ios_base::fmtflags flags = out.flags(ios::fixed);

    const double sigma = sqrt(m.variance());

    out << m.value() << " +/- " << sigma;

    out.precision(2);
    out << " ("  << 100 * sigma / m.value() << "%)";

    out.flags(flags);
    out.precision(precision);

    return out;
}

// Simple logic
//
Measurement s8::operator +(const Measurement &m1, const Measurement &m2)
{
    Measurement m = m1;

    m += m2;

    return m;
}

Measurement s8::operator -(const Measurement &m1, const Measurement &m2)
{
    Measurement m = m1;

    m -= m2;

    return m;
}

Measurement s8::operator *(const Measurement &m1, const Measurement &m2)
{
    Measurement m = m1;

    m *= m2;

    return m;
}

Measurement s8::operator /(const Measurement &m1, const Measurement &m2)
{
    Measurement m = m1;

    m /= m2;

    return m;
}

Measurement s8::operator %(const Measurement &m1, const Measurement &m2)
{
    Measurement m = m1;

    m %= m2;

    return m;
}

// Shortcuts: no temporary variable is created
//
void s8::operator+=(Measurement &m1, const Measurement &m2)
{
    m1.setValue(m1.value() + m2.value());
    m1.setVariance(m1.variance() + m2.variance());
}

void s8::operator-=(Measurement &m1, const Measurement &m2)
{
    m1.setValue(m1.value() - m2.value());
    m1.setVariance(m1.variance() + m2.variance());
}

void s8::operator*=(Measurement &m1, const Measurement &m2)
{
    m1.setValue(m1.value() * m2.value());
    m1.setVariance(pow(m2.value(), 2) * m1.variance() +
                   pow(m1.value(), 2) * m2.variance());
}

void s8::operator/=(Measurement &m1, const Measurement &m2)
{
    if (0 == m2.value())
        throw runtime_error("Failed to divide by null Measurement");

    m1.setValue(m1.value() / m2.value());
    m1.setVariance(pow(1. / m2.value(), 2) * m1.variance() +
                   pow(m1.value() / pow(m2.value(), 2), 2) *
                    m2.variance());
}

void s8::operator%=(Measurement &m1, const Measurement &m2)
{
    if (0 == m2.value())
        throw runtime_error("Failed to divide by null Measurement");

    const double ratio = m1.value() / m2.value();

    m1.setValue(ratio);
    m1.setVariance(((1 - 2 * ratio) * m1.variance() +
                   pow(ratio, 2) * m2.variance()) / pow(m2.value(), 2));
}
