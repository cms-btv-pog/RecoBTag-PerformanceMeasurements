
// PtrelSolver
// Author: Victor E. Bazterra, UIC (2008)

#ifndef PtrelByFitting_h
#define PtrelByFitting_h

#include "PtrelSolver.h"

class PtrelByFitting : public PtrelSolver
{

public:

    //! Default constructor.
    explicit PtrelByFitting(char const * filename) : PtrelSolver(filename) {}

    //! Measure efficiency
    virtual void solve(char const *, char const *);

private:

    //! Compute the efficiency
    bool compute(
        TH1 *,
        Flavor::Type dependency,
        ValueVector const &,
        CovarianceVector const &,
        ValueVector const &,
        CovarianceVector const &
    );

public:

    ClassDef(PtrelByFitting, 1)

};

#endif
