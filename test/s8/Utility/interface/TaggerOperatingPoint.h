/**
 * TaggerOperatingPoint
 * s8
 *
 * Created by Samvel Khalatian on Nov 15, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TAGGER_OPERATING_POINT
#define S8_TAGGER_OPERATING_POINT

#include <string>

#include "S8Tree/interface/S8Jet.h"

namespace s8
{
    class TaggerOperatingPoint
    {
        public:
            enum OperatingPoint {
                TCHET,  TCHEM,  TCHEL,
                TCHPT,  TCHPM,  TCHPL,
                SSVT,   SSVM,   SSVL,
                SSVHET, SSVHEM,
                SSVHPT
            };

            TaggerOperatingPoint() throw();
            ~TaggerOperatingPoint() throw();

            void initWithOperatingPoint(const OperatingPoint &);

            operator double() const;
            Jet::BTag btag() const;

        private:
            double    _operatingPoint;
            Jet::BTag _btag;

    };

    // Helpers
    //
    TaggerOperatingPoint &operator<<(TaggerOperatingPoint &,
                                     const std::string &);

    bool operator==(const TaggerOperatingPoint::OperatingPoint &,
                    const TaggerOperatingPoint &);

    bool operator==(const TaggerOperatingPoint &, const TaggerOperatingPoint &);
    bool operator!=(const TaggerOperatingPoint &, const TaggerOperatingPoint &);
}

#endif
