
// PtrelSolver
// Author: Victor E. Bazterra, UIC (2008)

#ifndef PtrelBySystem4_h
#define PtrelBySystem4_h

#include <vector>

#include "PtrelSolver.h"

class PtrelBySystem4 : public PtrelSolver
{

public:

    //! Default constructor.
    explicit PtrelBySystem4(char const * filename) : PtrelSolver(filename) {}

    //! Measure efficiency
    virtual void solve(char const *, char const *);
    
private:

    //! Histogram collection
    typedef std::vector<TH1D> HistogramCollection;

    //! Compute the coefficients between mc efficiencies
    bool makeCoefficients (
        TFile *,
        TFile *,
        TPRegexp &,
        TPRegexp &
    ) const;

    //! Compute the efficiency
    /* bool compute(
       TH1 *,
       ValueVector const &,
       TVectorD const &,
       TVectorD const &,
       ValueVector const &,
       CovarianceVector const &,
       ValueVector const &,
       CovarianceVector const &
    );*/

public:

    ClassDef(PtrelBySystem4, 1)

};

#endif
