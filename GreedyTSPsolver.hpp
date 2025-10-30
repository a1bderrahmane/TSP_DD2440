#include "ITSPsolver.hpp"
class GreedyTSPSolver : public ITSPSolver
{
public:
    vector<int> solve(vector<Point> &points) override;
};