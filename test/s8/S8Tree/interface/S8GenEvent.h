/**
 * GenEvent
 * s8 
 *
 * Created by Samvel Khalatian on Oct 11, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_GENEVENT
#define S8_GENEVENT

namespace s8
{
    class GenEvent
    {
        public:
            enum GluonSplitting { BB, CC, GLUON_SPLITTINGS};

            GenEvent() throw();

            GenEvent(const GenEvent &);
            GenEvent &operator =(const GenEvent &);

            void reset();

            bool isGluonSplitting() const;
            bool isGluonSplitting(const GluonSplitting &) const;
            double ptHat() const;



            void setGluonSplitting(const GluonSplitting &, const bool &);
            void setPtHat(const double &);

        private:
            // The same event may have several gluon splittings: bbbar, ccbar,
            // etc. Such events are supported even though they are highly
            // suppressed.
            //
            bool   _gluonSplittings[GLUON_SPLITTINGS];
            double _ptHat;
    };
}

#endif
