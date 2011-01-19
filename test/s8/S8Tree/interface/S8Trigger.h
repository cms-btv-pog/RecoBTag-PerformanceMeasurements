/**
 * Trigger
 * s8 
 *
 * Created by Samvel Khalatian on Oct 22, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TRIGGER
#define S8_TRIGGER

namespace s8
{
    class Trigger
    {
        public:
            typedef unsigned int Hash;

            Trigger() throw();

            Hash hash() const;
            int version() const;
            int prescale() const;

            virtual operator bool() const;

            void setHash(const Hash &);
            void setVersion(const int &);
            void setPrescale(const int &);
            void setIsPass(const bool &);

        private:
            Hash _hash;
            char _version;
            int  _prescale;
            bool _isPass;
    };
}

#endif
