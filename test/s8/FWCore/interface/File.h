/**
 * File, InputFile, OutputFile
 * core
 *
 * Created by Samvel Khalatian on Nov 12, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef CORE_FILE
#define CORE_FILE

#include <memory>
#include <string>

class TFile;

namespace core
{
    class File
    {
        public:
            enum Mode { READ = 0x1, RECREATE = 0x10 };

            File() throw();
            virtual ~File() throw();

            virtual std::string name() const;
            virtual void setName(const std::string &);

            virtual void open() = 0;
            virtual bool isOpen() const;
            virtual void close();

        protected:
            virtual void open(const Mode &);

            TFile *file() const;

        private:
            std::auto_ptr<TFile> _file;
            std::string          _name;
    };

    class InputFile: virtual public File
    {
        public:
            InputFile() throw();
            virtual ~InputFile() throw();

            virtual void open();
            virtual void process() = 0;
    };

    class OutputFile: virtual public File
    {
        public:
            OutputFile() throw();
            virtual ~OutputFile() throw();

            virtual void open();
    };
}

#endif
