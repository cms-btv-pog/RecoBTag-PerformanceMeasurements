/**
 * Jet
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_JET
#define S8_JET

#include <vector>

#include <TLorentzVector.h>

namespace s8
{
    class Jet
    {
        public:
            // Different BTaggers specification (to speed up search for
            // b-tagger). BTAGS is only for internal use. Any attempt to
            // use it outside will raise exception.
            //
            enum BTag { TCHE, TCHP, JP, SSV, SSVHE, SSVHP, BTAGS};

            Jet() throw();

            int flavour() const;
            int tracks() const;

            TLorentzVector &p4();
            const TLorentzVector &p4() const;

            double btag(const BTag &) const;



            void setFlavour(const int &);
            void setTracks(const int &);

            void setBTag(const BTag &, const double &);

        private:
            int _flavour;
            int _tracks;

            typedef std::vector<double> BTagCollection;

            TLorentzVector _p4;
            BTagCollection _btag;
    };

    inline int Jet::flavour() const
    {
        return _flavour;
    }

    inline int Jet::tracks() const
    {
        return _tracks;
    }

    inline TLorentzVector &Jet::p4()
    {
        return _p4;
    }

    inline const TLorentzVector &Jet::p4() const
    {
        return _p4;
    }

    inline double Jet::btag(const BTag &tag) const
    {
        // out_of_range exception will be raised if wrong BTag is specified.
        // BTAGS will be treated as ERROR
        //
        return _btag.at(tag);
    }



    inline void Jet::setFlavour(const int &flavour)
    {
        _flavour = flavour;
    }

    inline void Jet::setBTag(const BTag &tag, const double &value)
    {
        // out_of_range exception will be raised if wrong BTag is specified.
        // BTAGS will be treated as Error.
        //
        _btag.at(tag) = value;
    }
}

#endif
