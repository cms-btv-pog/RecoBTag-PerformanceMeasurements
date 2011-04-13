/**
 * PythiaOptionsDelegate
 * s8
 *
 * Created by Samvel Khalatian on Dec 16, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PYTHIA_OPTIONS_DELEGATE
#define S8_PYTHIA_OPTIONS_DELEGATE

namespace s8
{
    class Range;

    class PythiaOptionsDelegate
    {
        public:
            enum GluonSplitting { 
                // Do nothing
                //
                KEEP,

                // keep only specific gluon splitting events (KEEP = BB + CC)
                //
                BB, CC, ONLY,

                // remove specific gluon splitting events (REMOVE = BB + CC)
                //
                NO_BB, NO_CC, REMOVE,

                // enhance specifig gluon splitting events (ENHANCE = BB + CC)
                // [count every Nth event without gluon splitting]
                //
                ADD_BB, ADD_CC, ENHANCE
            };

            PythiaOptionsDelegate() throw();
            virtual ~PythiaOptionsDelegate() throw();

            virtual void optionGluonSplittingIsSet(const GluonSplitting &);
            virtual void optionPtHatIsSet(const Range &);
    };
}

#endif
