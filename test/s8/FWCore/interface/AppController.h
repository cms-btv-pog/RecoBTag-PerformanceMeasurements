/**
 * AppController
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_APPLICATION_CONTROLLER
#define CORE_APPLICATION_CONTROLLER

#include <string>
#include <vector>

#include "interface/GenericOptionsDelegate.h"

namespace core
{
    class Options;
    class InputFile;
    class OutputFile;

    // AppController inherits from the GenericOptionsDelegate for descendatns.
    //
    class AppController: public GenericOptionsDelegate
    {
        public:
            AppController() throw();
            virtual ~AppController() throw();

            // Parse Arguments
            //
            bool init(int argc, char *argv[]) throw();

            // Run the code
            //
            bool run() throw();

            // GenericOptionsDelegate interface
            //
            virtual void optionOutputIsSet(const std::string &);
            virtual void optionInputIsSet(const std::string &);

        protected:
            // Let descendants to Initialize. An exception should be raised
            // with description in case of any error.
            //
            virtual void init() = 0;

            virtual void applicationWillRun() = 0;
            virtual void applicationDidRun() = 0;

            virtual Options    *createOptions() = 0;
            virtual InputFile  *createInputFile() = 0;
            virtual OutputFile *createOutputFile() = 0;

        private:
            bool processInput() throw();

            typedef std::vector<std::string> InputFileNames;

            std::string    _outputFileName;
            InputFileNames _inputFileNames;
    };
}

#endif
