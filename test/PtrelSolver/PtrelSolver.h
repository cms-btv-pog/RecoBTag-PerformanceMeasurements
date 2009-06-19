
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
        useTwoTemplates();
        taggedLightTemplate_ = true;
    }

    //! Destructor.
    virtual ~PtrelSolver() {}

    //! Reset the solver
    void reset()
    {
        templates_ = 0;
        rebin_ = std::vector<Int_t>(Dependency::Dimension, 1);
    }

    //! Set the fitter for using three templates
    void useTwoTemplates()
    {
        Info(__FUNCTION__, "Using two templates for ptrel fitting.");
        fitFlavors_.clear();
        fitFlavors_.push_back(Flavor::b);
        fitFlavors_.push_back(Flavor::cl);
    }

    //! Set the fitter for using three templates
    void useThreeTemplates()
    {
        Info(__FUNCTION__, "Using three templates for ptrel fitting.");
        fitFlavors_.clear();
        fitFlavors_.push_back(Flavor::b);
        fitFlavors_.push_back(Flavor::c);
        fitFlavors_.push_back(Flavor::l);
    }

    //! Set the use of tagged light templates
    void setTaggedLightTemplate(bool flag = true)
    {
        taggedLightTemplate_ = flag;
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

    //! Process 2D histogram and rebinning if it is needed
    TH2 * processTH2(TObject *) const;

    static char const * directory;

    TFile * templates_;

private:

    Fit::Type fittype_;

    std::vector<Int_t> rebin_;

    bool taggedLightTemplate_;

    std::vector<Flavor::Type> fitFlavors_;

    //! Measure the flavor content for a sample
    bool measure(TFile *, TH2 *, ValueVector &, ValueVector &);

    //! Ptrel fit (by function) of a given sample
    bool fit(TH1 *, TF1 *, TVectorD &, TVectorD &);

    //! Ptrel fit (by histograms) of a given sample
    bool fit(TFile *, TH1 *, TObjArray *, TVectorD &, TVectorD &);

    TF1 * combinedFunctions_;

    TObjArray * combinedHistograms_;

    //! Return a combined function by given the sample name and bin number
    bool combinedTemplates(const char *, Int_t);

public:

    ClassDef(PtrelSolver, 1)
};

#endif

