/**
 * File, InputFile, OutputFile
 * s8
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_FILE
#define S8_FILE

#include "FWCore/interface/File.h"

namespace s8
{
    class InputFileDelegate;
    class OutputFileDelegate;

    class InputFile: public core::InputFile
    {
        public:
            InputFile() throw();
            virtual ~InputFile() throw();

            InputFileDelegate *delegate() const;
            void setDelegate(InputFileDelegate *);

            // core::InputFile interface
            //
            virtual void open();
            virtual void process();
            virtual void close();

        private:
            InputFileDelegate *_delegate;
    };

    class OutputFile: public core::OutputFile
    {
        public:
            OutputFile() throw();
            virtual ~OutputFile() throw();

            OutputFileDelegate *delegate() const;
            void setDelegate(OutputFileDelegate *);

            // core::OutputFile interface
            //
            virtual void open();
            virtual void close();

        private:
            OutputFileDelegate *_delegate;
    };
}

#endif
