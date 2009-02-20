
#ifndef PtrelSolver_h
#define PtrelSolver_h

#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPRegexp.h"
#include "TString.h"
#include "TVectorD.h"

#include "PtrelSolverDependencies.h"

class PtrelSolver : public TObject
{

public:

    //! Default constructor.
    explicit PtrelSolver(char const * filename, Fit::Type fittype)
    {
        reset();
        fittype_ = fittype;
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
    typedef std::vector<std::vector<TVectorD> > ValueMatrix;

    //! Measure the flavor content for a given sample
    bool measure(TFile *, TFile *, char const *, ValueVector &, ValueVector &);

    //! Measure the flavor content for samples matching a pattern
    bool measure(TFile *, TFile *, TPRegexp, StringVector &, ValueMatrix &, ValueMatrix &);

    //! Process 2D histogram and rebinning if it is needed
    TH2 * readTH2(TFile *, const char *) const;

    static char const * directory;

    TFile * templates_;

private:

    Fit::Type fittype_;

    std::vector<Int_t> rebin_;

    //! Measure the flavor content for a sample
    bool measure(TFile *, TH2 *, ValueVector &, ValueVector &);

    //! Process 2D histogram and rebinning if it is needed
    TH2 * processTH2(TObject *) const;

    //! Ptrel fit (by function) of a given sample
    bool fit(TH1 *, TF1 *, TVectorD &, TVectorD &);

    //! Ptrel fit (by histograms) of a given sample
    bool fit(TFile *, TH1 *, TObjArray *, TVectorD &, TVectorD &);

    //! Return a combined function by given the sample name and bin number
    TF1 * combinedFunction(const char *, Int_t);

    //! Return a set of histograms by given the sample name and bin number
    TObjArray * combinedHistograms(const char *, Int_t);

public:

    ClassDef(PtrelSolver, 1)
};

#endif

