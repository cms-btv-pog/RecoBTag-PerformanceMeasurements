#ifndef PtrelSolver_h
#define PtrelSolver_h

#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TObject.h"
#include "TPRegexp.h"
#include "TString.h"
#include "TVectorD.h"

#include "PtrelSolverDependencies.h"

class PtrelSolver : public TObject
{

public:

    //! Default constructor.
    explicit PtrelSolver(char const * filename)
    {
        reset();
        templates(filename);
    }

    //! Destructor.
    virtual ~PtrelSolver() {}

    //! Reset the solver
    void reset()
    {
        templates_ = 0;
        rebin_ = std::vector<Int_t>(Dependency::Dimension, 1);
    }

    //! Load templates
    void templates(char const * filename);

    //! Measure efficiency
    virtual void solve(char const *, char const *) const {}

protected:

    typedef std::vector<TString> StringVector;
    typedef std::vector<TVectorD> ValueVector;
    typedef std::vector<TMatrixD> CovarianceVector;

    typedef std::vector<std::vector<TVectorD> > ValueMatrix;
    typedef std::vector<std::vector<TMatrixD> > CovarianceMatrix;

    //! Measure the flavor content for a given sample
    bool measure(TFile *, TFile *, char const *, ValueVector &, CovarianceVector &);

    //! Measure the flavor content for samples matching a pattern
    bool measure(TFile *, TFile *, TPRegexp, StringVector &, ValueMatrix &, CovarianceMatrix &);

    //! Process 2D histogram and rebinning if it is needed
    TH2 * readTH2(TFile *, const char *) const;

    static char const * directory;

private:

    TFile * templates_;

    std::vector<Int_t> rebin_;

    //! Measure the flavor content for a samples
    bool measure(TFile *, TH2 *, ValueVector &, CovarianceVector &);

    //! Process 2D histogram and rebinning if it is needed
    TH2 * processTH2(TObject *) const;

    //! Ptrel fit of a given sample
    bool fit(TH1 *, TF1 *, TVectorD &, TMatrixD &);

    //! Return a combined function by given the sample name and bin number
    TF1 * combinedFunction(const char *, Int_t);

public:

    ClassDef(PtrelSolver, 1)
};

#endif

