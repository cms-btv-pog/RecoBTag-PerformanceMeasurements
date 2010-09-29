/**
 * Tools
 * s8
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_TOOLS
#define S8_TOOLS

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

class TLorentzVector;
class TVector3;

namespace s8
{
    namespace tools
    {
        void setP4(TLorentzVector &,
                   const math::XYZTLorentzVector &);

        void setVertex(TVector3 &,
                       const math::XYZPoint &);
    }
}

#endif
