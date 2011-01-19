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
    class PythiaOptionsDelegate
    {
        public:
            enum GluonSplitting { KEEP,     // keep all events including the
                                            // Gluon Splitting

                                  REMOVE,   // remove events with gluon splitting

                                  ONLY,     // keep only events with the
                                            // Gluon Splitting

                                  ENHANCE,  // Enhanse sample with Gluon
                                            // Splitting: count the non-gluon
                                            // splitting event events every Nth

                                  SIMULATE  // Simulate Gluon splitting: 2 jets
                                            // with high pt and deltaR
            };

            PythiaOptionsDelegate() throw();
            virtual ~PythiaOptionsDelegate() throw();

            virtual void optionGluonSplittingIsSet(const GluonSplitting &);
            virtual void optionPtHatIsSet(const int &min, const int &max);
    };
}

#endif
