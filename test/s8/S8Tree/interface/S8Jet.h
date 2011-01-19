/**
 * Jet
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_JET
#define S8_JET

class TLorentzVector;

namespace s8
{
    class GenParticle;

    class Jet
    {
        public:
            // Different BTaggers specification (to speed up search for
            // b-tagger). BTAGS are only for internal use.
            //
            enum BTag { TCHE, TCHP, JP, SSV, SSVHE, SSVHP, BTAGS};

            Jet() throw();
            ~Jet() throw();

            Jet(const Jet &);
            Jet &operator =(const Jet &);

            void reset();

            int flavour() const;
            int tracks() const;

            TLorentzVector *p4();
            const TLorentzVector *p4() const;

            GenParticle *genParticle();
            const GenParticle *genParticle() const;

            double btag(const BTag &) const;



            void setFlavour(const int &);
            void setTracks(const int &);

            void setBTag(const BTag &, const double &);

        private:
            int _flavour;
            int _tracks;

            TLorentzVector *_p4;
            GenParticle    *_genParticle;

            // Low level array is used instead of vector to speed up
            // compilation
            //
            double _btag[BTAGS];
    };
}

#endif
