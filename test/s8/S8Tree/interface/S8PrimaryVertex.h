/**
 * PrimaryVertex
 * s8 
 *
 * Created by Samvel Khalatian on Sep 28, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_PRIMARY_VERTEX
#define S8_PRIMARY_VERTEX

#include <memory>

class TVector3;

namespace s8
{
    class PrimaryVertex
    {
        public:
            PrimaryVertex() throw();
            ~PrimaryVertex() throw();

            PrimaryVertex(const PrimaryVertex &);
            PrimaryVertex &operator =(const PrimaryVertex &);

            void reset();

            TVector3 *vertex();
            const TVector3 *vertex() const;

            double ndof() const;
            double rho() const;



            void setNdof(const double &);
            void setRho(const double &);

        private:
            TVector3 *_vertex;

            double   _ndof;
            double   _rho;
    };
}

#endif
