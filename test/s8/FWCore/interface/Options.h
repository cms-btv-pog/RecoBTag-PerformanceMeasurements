/**
 * Options
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_OPTIONS
#define CORE_OPTIONS

#include <iosfwd>

namespace boost
{
    namespace program_options
    {
        class options_description;
    }
}

namespace po = boost::program_options;

namespace core
{
    class Options
    {
        public:
            Options() throw();
            virtual ~Options() throw();

            // Initialize object: construction options, etc.
            //
            virtual void init() = 0;

            // Access created options_description
            //
            virtual po::options_description *description() const = 0;

            // Print set options in output stream
            //
            virtual void print(std::ostream &) const = 0;
    };

    std::ostream &operator<<(std::ostream &, const Options &);
}

#endif
