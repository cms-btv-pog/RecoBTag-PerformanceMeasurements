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
            ~EventID() throw();

            void reset();

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
}

#endif
