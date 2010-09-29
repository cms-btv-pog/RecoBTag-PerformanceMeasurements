/**
 * GenParticle
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_GENPARTICLE
#define S8_GENPARTICLE

#include <TLorentzVector.h>
#include <TVector3.h>

namespace s8
{
    class GenParticle
    {
        public:
            GenParticle() throw();

            int id() const;
            int parentId() const;

            TLorentzVector &p4();
            const TLorentzVector &p4() const;

            TVector3 &vertex();
            const TVector3 &vertex() const;



            void setId(const int &);
            void setParentId(const int &);

        private:
            int _id;
            int _parentId;

            TLorentzVector  _p4;
            TVector3        _vertex;
    };

    inline int GenParticle::id() const
    {
        return _id;
    }

    inline int GenParticle::parentId() const
    {
        return _parentId;
    }

    inline TLorentzVector &GenParticle::p4()
    {
        return _p4;
    }

    inline const TLorentzVector &GenParticle::p4() const
    {
        return _p4;
    }

    inline TVector3 &GenParticle::vertex()
    {
        return _vertex;
    }

    inline const TVector3 &GenParticle::vertex() const
    {
        return _vertex;
    }



    inline void GenParticle::setId(const int &id)
    {
        _id = id;
    }

    inline void GenParticle::setParentId(const int &id)
    {
        _parentId = id;
    }
}

#endif
