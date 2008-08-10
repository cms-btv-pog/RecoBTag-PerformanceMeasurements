#ifndef PtrelSolverDependencies_h
#define PtrelSolverDependencies_h

#include "TObject.h"

//! Define the set of dependencies
struct Dependency : public TObject
{
    enum Type { ptrel = 0, pT, eta};
    static const char * Label[];
    static const char * Name[];
    static const long Dimension;
    ClassDef(Dependency, 1)
};

//! Define the set of dependencies
struct Flavor : public TObject
{
    enum Type { b = 0, cl};
    static const char * Label[];
    static const char * Name[];
    static const long Dimension;
    ClassDef(Flavor, 1)
};

#endif
