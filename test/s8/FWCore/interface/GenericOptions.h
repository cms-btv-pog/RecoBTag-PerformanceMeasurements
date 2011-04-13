/**
 * GenericOptions
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_GENERIC_OPTIONS
#define CORE_GENERIC_OPTIONS

#include <iosfwd>
#include <string>
#include <vector>
#include <utility>

#include "interface/Options.h"

namespace boost
{
    namespace program_options
    {
        class positional_options_description;
    }
}

namespace core
{
    class GenericOptionsDelegate;

    class GenericOptions: public Options
    {
        public:
            GenericOptions() throw();
            virtual ~GenericOptions() throw();

            virtual void init();

            GenericOptionsDelegate *delegate() const;
            void setDelegate(GenericOptionsDelegate *);

            virtual po::options_description *description() const;
            po::positional_options_description *positional() const;

            virtual void print(std::ostream &) const;

        private:
            void optionDebugIsSet(const std::string &);
            void optionOutputIsSet(const std::string &);
            void optionEventsIsSet(const int &);
            void optionSkipEventsIsSet(const int &);
            void optionInputIsSet(const std::vector<std::string> &);

            void readInputFilesFromTXT(const std::string &);

            GenericOptionsDelegate                            *_delegate;
            std::auto_ptr<po::options_description>             _description;
            std::auto_ptr<po::positional_options_description>  _positional;

            std::string              _debug;
            std::string              _output;
            int                      _events;
            int                      _skip_events;

            // pirs: (Filename,DoesFileExist)
            //
            typedef std::pair<std::string, bool> FilePair;
            typedef std::vector<FilePair> Files;
            Files _inputs;
    };
}

#endif
