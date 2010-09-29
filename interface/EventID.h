/**
 * EventID
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_EVENTID
#define S8_EVENTID

namespace s8
{
    class EventID
    {
        public:
            EventID() throw();

            int run() const;
            int lumiBlock() const;
            int event() const;



            void setRun(const int &);
            void setLumiBlock(const int &);
            void setEvent(const int &);

        private:
            int _run;
            int _lumiBlock;
            int _event;
    };

    inline int EventID::run() const
    {
        return _run;
    }

    inline int EventID::lumiBlock() const
    {
        return _lumiBlock;
    }

    inline int EventID::event() const
    {
        return _event;
    }
}

#endif
