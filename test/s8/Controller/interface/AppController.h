/**
 * AppController
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_APPLICATION_CONTROLLER
#define S8_APPLICATION_CONTROLLER

#include <ctime>
#include <memory>

#include "FWCore/interface/AppController.h"
#include "IO/interface/FileDelegate.h"

namespace core
{
    class Debug;
}

namespace s8
{
    class Analyzer;
    class InputFile;
    class OutputFile;

    class AppController: public core::AppController,
                         public InputFileDelegate,
                         public OutputFileDelegate
    {
        public:
            AppController() throw();
            virtual ~AppController() throw();

        protected:
            // core::AppController interface
            //
            virtual void init();

            virtual void applicationWillRun();
            virtual void applicationDidRun();

            // descendants interface
            //
            virtual Analyzer *createAnalyzer() = 0;

            // Generic Options Delegate interface
            //
            virtual void optionDebugIsSet(const std::string &);
            virtual void optionEventsIsSet(const int &);

        protected:
            // core::AppController interface
            //
            virtual core::Options    *createOptions();
            virtual core::InputFile  *createInputFile();
            virtual core::OutputFile *createOutputFile();

            // InputFileDelegate interface
            //
            virtual bool inputFileShouldOpen(const std::string &);
            virtual void inputFileDidOpen(TFile *);
            virtual bool inputFileShouldLoadJets();
            virtual bool inputFileShouldLoadMuons();
            virtual bool inputFileShouldLoadPrimaryVertices();
            virtual void inputFileDidLoadEvent(const Event *);
            virtual bool inputFileShouldContinue();
            virtual void inputFileWillClose(TFile *);

            // OutputFileDelegate interface
            //
            virtual void outputFileDidOpen(TFile *);
            virtual void outputFileWillClose(TFile *);

        private:
            Analyzer *_analyzer;
            std::auto_ptr<InputFile>  _inputFile;
            std::auto_ptr<OutputFile> _outputFile;

            std::auto_ptr<core::Debug> _debug;

            int _maxEvents;
            int _processedEvents;

            clock_t _startClock;
    };
}

#endif
