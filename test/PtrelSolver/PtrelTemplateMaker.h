
#ifndef PtrelTemplateMaker_h
#define PtrelTemplateMaker_h

#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TPRegexp.h"
#include "TObject.h"
#include "TString.h"

#include "PtrelUtils.h"
#include "PtrelSolverDependencies.h"

class PtrelTemplateMaker : public TObject
{

public:

    //! Constructor by function.
    explicit PtrelTemplateMaker()
    {
        reset();
    }

    //! Destructor.
    virtual ~PtrelTemplateMaker() {}

    //! Reset to the default configuration.
    void reset()
    {
        rebin_ = std::vector<Int_t>(Dependency::Dimension, 1);
        functions_ = std::vector<TF1>(Flavor::Dimension, TF1());
    }

    //! Set rebinning options
    void function(Dependency::Type dependency, TF1 const & form)
    {
        functions_[dependency] = form;
    }

    //! Set rebinning options
    void rebin(Dependency::Type dependency, Int_t value)
    {
        rebin_[dependency] = value;
    }

    //! Make all the templates
    void make(char const *, char const *) const;

private:

    std::vector<TF1> functions_;

    std::vector<Int_t> rebin_;

    static char const * directory;

    TH2 * processTH2(TObject *) const;

    //! Implementation for creating templates
    bool makeTemplates(TFile *, TFile *) const;

    //! Implementation for creating efficiencies
    bool makeEfficiencies (TFile *, TFile *, TPRegexp, TPRegexp) const;

    //! Implementation for creating coefficiencies between efficiencies
    bool makeCoefficients (TFile *, TPRegexp, TPRegexp) const;

public:

    ClassDef(PtrelTemplateMaker, 1)

};

#endif
