#include "ITSPsolver.hpp"
class GreedyTSPSolver : public ITSPSolver
{
public:
    vector<int> solve(const vector<Point> &points) override;
};