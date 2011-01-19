/**
 * Lepton
 * s8 
 *
 * Created by Samvel Khalatian on Oct 20, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_LEPTON
#define S8_LEPTON

#include <utility>

class TLorentzVector;
class TVector3;

namespace s8
{
    class GenParticle;

    class Lepton
    {
        public:
            // Value and associated Error pair
            //
            typedef std::pair<double, double> ImpactParameter;

            Lepton() throw();
            ~Lepton() throw();

            // Take care of copying b/c of pointers
            //
            Lepton(const Lepton &);
            Lepton &operator =(const Lepton &);

            void reset();

            ImpactParameter *impactParameter();
            const ImpactParameter *impactParameter() const;

            TLorentzVector *p4();
            const TLorentzVector *p4() const;

            TVector3 *vertex();
            const TVector3 *vertex() const;

            GenParticle *genParticle();
            const GenParticle *genParticle() const;

        private:
            ImpactParameter *_impactParameter;
            TLorentzVector  *_p4;
            TVector3        *_vertex;
            GenParticle     *_genParticle;
    };
}

#endif
