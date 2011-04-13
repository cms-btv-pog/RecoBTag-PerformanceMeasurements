/**
 * Range
 * s8
 *
 * Created by Samvel Khalatian on Feb 03, 2011
 * Copyright 2010, All rights reserved
 */

#ifndef S8_RANGE
#define S8_RANGE

#include <string>
#include <iosfwd>

namespace s8
{
    /*
     * Unchecked Range: no check is done if MIN < MAX
     */
    class Range
    {
        public:
            Range() throw();
            ~Range() throw();

            double minimum() const throw();
            double maximum() const throw();

            void setMinimum(const double &);
            void setMaximum(const double &);

        private:
            double _minimum;
            double _maximum;
    };

    // Helpers
    // -----------------------------------------------------------------------

    // Test if value is in range: [min, max)
    // Note: range is tested for validity
    //
    bool isValueInRange(const double &, const Range &);

    // Parse string of format: MIN..MAX
    // (doubles are accepted)
    //
    void parse(Range &, const std::string &);

    std::ostream &operator<<(std::ostream &, const Range &);
}

#endif
