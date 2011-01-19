/**
 * InputFileDelegate
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_FILE_DELEGATE
#define S8_FILE_DELEGATE

#include <string>

class TFile;

namespace s8
{
    class Event;

    class InputFileDelegate
    {
        public:
            virtual ~InputFileDelegate() throw();

            virtual bool inputFileShouldOpen(const std::string &);
            virtual void inputFileDidOpen(TFile *);
            virtual void inputFileDidLoadEvent(const Event *);
            virtual bool inputFileShouldLoadJets();
            virtual bool inputFileShouldLoadElectrons();
            virtual bool inputFileShouldLoadMuons();
            virtual bool inputFileShouldLoadPrimaryVertices();
            virtual bool inputFileShouldLoadTriggers();
            virtual bool inputFileShouldContinue();
            virtual void inputFileWillClose(TFile *);
    };

    class OutputFileDelegate
    {
        public:
            virtual ~OutputFileDelegate() throw();

            virtual void outputFileDidOpen(TFile *);
            virtual void outputFileWillClose(TFile *);
    };
}

#endif
