#ifndef S8_SOLUTION
#define S8_SOLUTION

typedef std::map<std::string, Measurement> Solution;

struct SolutionInBin
{
    Solution     solution;
    Measurement  bin;
    unsigned int binID;
};

#endif
