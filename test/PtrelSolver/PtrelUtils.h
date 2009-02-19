
#ifndef PtrelUtils_h
#define PtrelUtils_h

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

#endif
