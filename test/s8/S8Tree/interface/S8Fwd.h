/**
 * Forward declarations
 * s8 
 *
 * Created by Samvel Khalatian on Nov 20, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_FORWARD_DECLARATIONS
#define S8_FORWARD_DECLARATIONS

#include <vector>

namespace s8
{
    class Jet;
    class Lepton;
    class PrimaryVertex;
    class Trigger;

    typedef std::vector<Jet *>           Jets;
    typedef std::vector<Lepton *>        Leptons;
    typedef std::vector<PrimaryVertex *> PrimaryVertices;
    typedef std::vector<Trigger *>       Triggers;
}

#endif
