
#ifndef PtrelUtils_h
#define PtrelUtils_h

#include "TPRegexp.h"

#define CallSafely(code) \
if (!code) \
{ \
   Error(__FUNCTION__, "Function call %s fail", #code); \
   return; \
}

#define CallSafelyZero(code) \
if (!code) \
{ \
   Error(__FUNCTION__, "Function call %s fail", #code); \
   return 0; \
}

#define GetSafely(pointer, code) \
pointer = code; \
if (!pointer) \
{ \
   Error(__FUNCTION__, "Action %s fail", #code); \
   return; \
}


#define GetSafelyZero(pointer, code) \
pointer = code; \
if (!pointer) \
{ \
   Error(__FUNCTION__, "Action %s fail", #code); \
   return 0; \
}


#define CreateSafely(type, pointer, code) \
type * pointer = (type *) code; \
if (!pointer) \
{ \
   Error(__FUNCTION__, "Action %s fail", #code); \
   return; \
}


#define CreateSafelyZero(type, pointer, code) \
type * pointer = (type *) code; \
if (!pointer) \
{ \
   Error(__FUNCTION__, "Action %s fail", #code); \
   return 0; \
}

#include "TH1.h"

void ptrelHistogramSetup(TH1*);

void efficiencyHistogramSetup(TH1*);

bool containsIdentifier(const char* objName_, const char* id_);

#include "TVectorD.h"
#include "TMatrixD.h"

class ErrorPropagator
{
public:

    ErrorPropagator() : removeCorrelations_(false) {}

    ErrorPropagator(const TMatrixD & matrix) : removeCorrelations_(false)
    {
        transformation(matrix);
    }

    ErrorPropagator(Int_t size, const Double_t * matrix) : removeCorrelations_(false)
    {
        transformation( TMatrixD(size, size, matrix) );
    }

    void transformation(TMatrixD const &);

    void transformation(Int_t size, const Double_t * matrix)
    {
        transformation( TMatrixD(size, size, matrix) );
    }

    TVectorD const & correctedErrors(TVectorD const &);

    Double_t operator()(TVectorD const &, TVectorD const &);

    void removeCorrelations(bool flag = true)
    {
        removeCorrelations_ = flag;
    }

    TMatrixD const & matrix()
    {
        return matrix_;
    }

    TMatrixD const & fadamard()
    {
        return hadamard_;
    }

    TMatrixD const & invHadamard()
    {
        return inverseHadamard_;
    }

private:

    bool removeCorrelations_;
    TVectorD correctedErrors_;
    TMatrixD matrix_, hadamard_, inverseHadamard_;
};

#endif
