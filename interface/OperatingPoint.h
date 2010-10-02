/**
 * Operating Point
 * 
 *
 * Created by Samvel Khalatian on Oct 2, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_OPERATING_POINT
#define S8_OPERATING_POINT

#include <string>

namespace s8
{
    class OperatingPoint
    {
        public:
            enum OP { TCHET,  TCHEM,  TCHEL,
                      TCHPT,  TCHPM,  TCHPL,
                      SSVT,   SSVM,   SSVL,
                      SSVHET, SSVHEM,
                      SSVHPT };

            // TCHEM is the default operating point
            //
            OperatingPoint();

            void set(const OP &);

            operator double() const;

        private:
            double _value;
    };

    OperatingPoint &operator<<(OperatingPoint &, const std::string &);

    inline OperatingPoint::operator double() const
    {
        return _value;
    }
}

#endif
