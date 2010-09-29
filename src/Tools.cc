/**
 * Tools
 * s8
 *
 * Created by Samvel Khalatian on Sep 29, 2010
 * Copyright 2010, All rights reserved
 */

#include <TLorentzVector.h>
#include <TVector3.h>

#include "RecoBTag/PerformanceMeasurements/interface/Tools.h"

namespace s8tools = s8::tools;

void s8tools::setP4(TLorentzVector &p4,
                       const math::XYZTLorentzVector &cmsswP4)
{
    p4.SetPxPyPzE(cmsswP4.px(),
                  cmsswP4.py(),
                  cmsswP4.pz(),
                  cmsswP4.energy());
}

void s8tools::setVertex(TVector3 &vector, const math::XYZPoint &point)
{
    vector.SetXYZ(point.x(), point.y(), point.z());
}
