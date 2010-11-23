/**
 * S8NumericInput
 * 
 *
 * Created by Samvel Khalatian on Nov 10, 2010
 * Copyright 2010, All rights reserved
 */

#include <cmath>

#include <iomanip>

#include "S8NumericInput.h"

using std::endl;
using std::fixed;
using std::left;
using std::make_pair;
using std::setprecision;
using std::setw;

std::ostream &operator<<(std::ostream &out,
                         const Measurement &measurement)
{
    const double sigma = sqrt(measurement.second);

    return out << setprecision(4) << fixed << measurement.first << " +/- " << sigma
        << " ("  << 100 * sigma / measurement.first << "%)";
}

Measurement &operator+=(Measurement &in, const Measurement &measurement)
{
    in.first += measurement.first;
    in.second += measurement.second;

    return in;
}

Measurement operator +(const Measurement &m1,
                       const Measurement &m2)
{
    return make_pair(m1.first + m2.first,
                     m1.second + m2.second);
}

Measurement operator *(const double &scaleFactor,
                       const Measurement &measurement)
{
    return make_pair(scaleFactor * measurement.first,
                     pow(scaleFactor, 2) * measurement.second);
}

Measurement operator /(const Measurement &m1,
                       const Measurement &m2)
{
    return make_pair(m1.first / m2.first,
                     pow(1. / m2.first, 2) * m1.second +
                        pow(m1.first / pow(m2.first, 2), 2) * m2.second);
}


Measurement operator %(const Measurement &m1,
                       const Measurement &m2)
{
    const double eff = m1.first / m2.first;

    return make_pair(eff,
                     ((1 - 2 * eff) * m1.second +
                        pow(eff, 2) * m2.second) / pow(m2.first, 2));
}

std::ostream &operator<<(std::ostream &out, const NumericInput &input)
{
    out << " [+] " << setw(15) << left << "n"
        << input.n.all << endl;
    out << " [+] " << setw(15) << left << "n mu"
        << input.n.mu << endl;
    out << " [+] " << setw(15) << left << "n tag"
        << input.n.tag << endl;
    out << " [+] " << setw(15) << left << "n mu tag"
        << input.n.muTag << endl;
    out << endl;

    out << " [+] " << setw(15) << left << "p"
        << input.p.all << endl;
    out << " [+] " << setw(15) << left << "p mu"
        << input.p.mu << endl;
    out << " [+] " << setw(15) << left << "p tag"
        << input.p.tag << endl;
    out << " [+] " << setw(15) << left << "p mu tag"
        << input.p.muTag << endl;

    return out;
}

std::ostream &operator<<(std::ostream &out, const FlavouredNumericInput &input)
{
    const numeric::FlavouredInputGroup &n = input.n;
    const numeric::FlavouredInputGroup &p = input.p;

    out << " [+] " << setw(15) << left << "n b" << n.b << endl;
    out << " [+] " << setw(15) << left << "n cl" << n.cl << endl;

    out << " [+] " << setw(15) << left << "n tag b" << n.tag.b << endl;
    out << " [+] " << setw(15) << left << "n tag cl" << n.tag.cl << endl;

    out << " [+] " << setw(15) << left << "n mu b" << n.mu.b << endl;
    out << " [+] " << setw(15) << left << "n mu cl" << n.mu.cl << endl;

    out << " [+] " << setw(15) << left << "n muTag b" << n.muTag.b << endl;
    out << " [+] " << setw(15) << left << "n muTag cl" << n.muTag.cl << endl;
    out << endl;

    out << " [+] " << setw(15) << left << "p b" << p.b << endl;
    out << " [+] " << setw(15) << left << "p cl" << p.cl << endl;

    out << " [+] " << setw(15) << left << "p tag b" << p.tag.b << endl;
    out << " [+] " << setw(15) << left << "p tag cl" << p.tag.cl << endl;

    out << " [+] " << setw(15) << left << "p mu b" << p.mu.b << endl;
    out << " [+] " << setw(15) << left << "p mu cl" << p.mu.cl << endl;

    out << " [+] " << setw(15) << left << "p muTag b" << p.muTag.b << endl;
    out << " [+] " << setw(15) << left << "p muTag cl" << p.muTag.cl << endl;
    out << endl;

    return out;
}

std::ostream &operator<<(std::ostream &out, const FlavouredEffGroup &group)
{
    out << " [+] " << setw(15) << left << "eff_tag_b" << group.tag.b << endl;
    out << " [+] " << setw(15) << left << "eff_tag_cl" << group.tag.cl << endl;
    out << " [+] " << setw(15) << left << "eff_mu_b" << group.mu.b << endl;
    out << " [+] " << setw(15) << left << "eff_mu_cl" << group.mu.cl << endl;

    return out;
}

std::ostream &operator<<(std::ostream &out, const Coefficients &coefficients)
{
    out << " [+] " << setw(15) << left << "alpha"
        << coefficients.alpha << endl;
    out << " [+] " << setw(15) << left << "beta"
        << coefficients.beta << endl;
    out << " [+] " << setw(15) << left << "gamma"
        << coefficients.gamma << endl;
    out << " [+] " << setw(15) << left << "delta"
        << coefficients.delta << endl;
    out << " [+] " << setw(15) << left << "kappaCL"
        << coefficients.kappaCL << endl;
    out << " [+] " << setw(15) << left << "kappaB"
        << coefficients.kappaB << endl;

    return out;
}

std::ostream &operator<<(std::ostream &out, const NumericInputGroup &group)
{
    out << " Bin" << endl;
    out << "     " << setw(15) << left << " "
        << group.bin << endl;
    out << endl;

    out << " Input" << endl;
    out << group.input << endl;

    out << " Flavoured Input" << endl;
    out << group.flavouredInput << endl;

    out << " Efficiencies" << endl;
    out << group.efficiency << endl;

    out << " Coefficients" << endl;
    out << group.coefficients << endl;

    return out;
}
