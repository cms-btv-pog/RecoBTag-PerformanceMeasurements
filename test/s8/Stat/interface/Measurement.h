/**
 * Measurement
 * s8
 *
 * Created by Samvel Khalatian on Dec 20, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_MEASUREMENT
#define S8_MEASUREMENT

#include <iosfwd>

namespace s8
{
    class Measurement
    {
        public:
            Measurement() throw();

            // Measurement(central_value, variance)
            //
            explicit Measurement(const double &, const double &) throw();

            double value() const;
            double variance() const;

            void setValue(const double &);
            void setVariance(const double &);

        private:
            double _value;
            double _variance;
    };

    std::ostream &operator <<(std::ostream &, const Measurement &);

    // Define simple logic
    //
    // Note: Measurements are assumed to be indenepdent in the error
    //       propagation
    //
    Measurement operator +(const Measurement &, const Measurement &);
    Measurement operator -(const Measurement &, const Measurement &);
    Measurement operator *(const Measurement &, const Measurement &);

    // Regular error propagation
    //
    Measurement operator /(const Measurement &, const Measurement &);

    // Binomial error propagation
    //
    Measurement operator %(const Measurement &, const Measurement &);

    // Shortcuts for the above operators
    //
    void operator+=(Measurement &, const Measurement &);
    void operator-=(Measurement &, const Measurement &);
    void operator*=(Measurement &, const Measurement &);

    // Regular error propagation
    //
    void operator/=(Measurement &, const Measurement &);

    // Binomial error propagation
    //
    void operator%=(Measurement &, const Measurement &);
}

#endif
