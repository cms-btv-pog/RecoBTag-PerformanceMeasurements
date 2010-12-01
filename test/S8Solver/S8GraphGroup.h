/**
 * S8GraphGroup
 * 
 *
 * Created by Samvel Khalatian on Nov 23, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_GRAPH_GROUP
#define S8_GRAPH_GROUP

#include <map>
#include <memory>
#include <stack>
#include <vector>

#include "S8NumericInput.h"
#include "S8Solution.h"

class TGraph;
class TGraphErrors;
class TObject;

typedef std::vector<NumericInputGroup> BinnedNumericInputGroup;
typedef std::vector<SolutionInBin> BinnedSolution;

struct Graph
{
    enum Type { PT, ETA, PHI };
};

struct FlavouredEffGraphGroup
{
    explicit FlavouredEffGraphGroup(const int &size);
    FlavouredEffGraphGroup(const BinnedSolution &);
    FlavouredEffGraphGroup(const BinnedNumericInputGroup &);

    void init(const int &, const bool &isMonteCarlo = false);

    std::auto_ptr<TGraphErrors> b;
    std::auto_ptr<TGraphErrors> cl;
};

struct EffGraphGroup
{
    explicit EffGraphGroup(const Graph::Type &, const int &);
    EffGraphGroup(const Graph::Type &, const BinnedSolution &);
    EffGraphGroup(const Graph::Type &, const BinnedNumericInputGroup &);

    void save(TDirectory *);

    FlavouredEffGraphGroup mu;
    FlavouredEffGraphGroup tag;
};

struct EffGraph
{
    EffGraph(const Graph::Type &,
             const BinnedNumericInputGroup &, const BinnedSolution &);
    ~EffGraph();

    void draw();
    void save(TDirectory *);

    EffGraphGroup mc;
    EffGraphGroup s8;
    EffGraphGroup scale;

    private:
        typedef std::stack<TObject *> Heap;

        Heap        _heaps;
        Graph::Type _type;
};

struct InputGraph
{
    InputGraph(const Graph::Type &, const int &size);

    std::auto_ptr<TGraphErrors> all;
    std::auto_ptr<TGraphErrors> mu;
    std::auto_ptr<TGraphErrors> tag;
    std::auto_ptr<TGraphErrors> muTag;
};

struct InputGraphGroup
{
    InputGraphGroup(const Graph::Type &, const BinnedNumericInputGroup &);
    ~InputGraphGroup();

    void draw();
    void save(TDirectory *);

    InputGraph n;
    InputGraph p;

    private:
        typedef std::stack<TObject *> Heap;

        Heap        _heaps;
        Graph::Type _type;
};

struct GraphGroup
{
    GraphGroup(const BinnedNumericInputGroup &, const BinnedSolution &,
               const Graph::Type & = Graph::PT);
    ~GraphGroup();

    void save(TDirectory *);
    void draw();

    // Note: objects will be automatically destroyed. Clone graph errors
    // if TMultiGraph is used
    //
    std::auto_ptr<TGraphErrors> alpha;
    std::auto_ptr<TGraphErrors> beta;
    std::auto_ptr<TGraphErrors> gamma;
    std::auto_ptr<TGraphErrors> delta;
    std::auto_ptr<TGraphErrors> kappaB;
    std::auto_ptr<TGraphErrors> kappaCL;

    EffGraph efficiency;

    InputGraphGroup input;

    private:
        typedef std::stack<TObject *> Heap;

        Heap _heaps;
        Graph::Type _type;
};

#endif
