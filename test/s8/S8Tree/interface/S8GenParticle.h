/**
 * GenParticle
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_GENPARTICLE
#define S8_GENPARTICLE

class TLorentzVector;
class TVector3;

namespace s8
{
    class GenParticle
    {
        public:
            GenParticle() throw();
            ~GenParticle() throw();

            // Take care of copying b/c of pointers
            //
            GenParticle(const GenParticle &);
            GenParticle &operator =(const GenParticle &);

            void reset();

            int id() const;
            int parentId() const;
            int status() const;

            TLorentzVector *p4();
            const TLorentzVector *p4() const;

            TVector3 *vertex();
            const TVector3 *vertex() const;



            void setId(const int &);
            void setParentId(const int &);
            void setStatus(const int &);

        private:
            int _id;
            int _parentId;
            int _status;

            TLorentzVector *_p4;
            TVector3       *_vertex;
    };
}

#endif
