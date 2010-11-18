/**
 * S8SolverInput
 * 
 *
 * Created by Samvel Khalatian on Nov 11, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_SOLVERINPUT
#define S8_SOLVERINPUT

class TH1;

namespace solver
{
    struct PlotGroup
    {
        PlotGroup() throw();

        TH1 *all;
        TH1 *tag;
        TH1 *mu;
        TH1 *muTag;
    };



    struct FlavouredPlot
    {
        FlavouredPlot() throw();
        virtual ~FlavouredPlot() throw();

        TH1 *b;
        TH1 *cl;
    };

    struct FlavouredPlotGroup: public FlavouredPlot
    {
        virtual ~FlavouredPlotGroup() throw();

        FlavouredPlot tag;
        FlavouredPlot mu;
        FlavouredPlot muTag;
    };
}

struct SolverInput
{
    solver::PlotGroup n;
    solver::PlotGroup p;
};

struct FlavouredSolverInput
{
    solver::FlavouredPlotGroup n;
    solver::FlavouredPlotGroup p;
};

#endif
