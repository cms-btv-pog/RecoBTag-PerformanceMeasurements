#ifndef PtrelByCounting_h
#define PtrelByCounting_h

#include "PtrelSolver.h"

class PtrelByCounting : public PtrelSolver
{

public:

    //! Default constructor.
    explicit PtrelByCounting(char const * filename) : PtrelSolver(filename) {}

    //! Measure efficiency
    virtual void solve(char const *, char const *)
    {
        Info(__FUNCTION__, "This function in not valid use instead solve(input, mistag, output)");
    }

    //! Measure efficiency
    void solve(char const *, char const *, char const *);

private:

    //! Compute the efficiency
    bool compute(
        TH1 *,
        TH1 *,
        TVectorD const &,
        ValueVector const &,
        CovarianceVector const &,
        ValueVector const &,
        CovarianceVector const &
    );

public:

    ClassDef(PtrelByCounting, 1)

};

#endif
