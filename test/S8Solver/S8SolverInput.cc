/**
 * S8SolverInput
 * 
 *
 * Created by Samvel Khalatian on Nov 11, 2010
 * Copyright 2010, All rights reserved
 */

#include "S8SolverInput.h"

using solver::PlotGroup;
using solver::FlavouredPlot;
using solver::FlavouredPlotGroup;

PlotGroup::PlotGroup() throw()
{
    all = 0;
    tag = 0;
    mu = 0;
    muTag = 0;
}

FlavouredPlot::FlavouredPlot() throw()
{
    b = 0;
    cl = 0;
}

FlavouredPlot::~FlavouredPlot() throw()
{
}

FlavouredPlotGroup::~FlavouredPlotGroup() throw()
{
}
