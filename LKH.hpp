#include "ITSPsolver.hpp"
class LKH : public ITSPSolver
{
public:
    vector<int> solve(const vector<Point> &points) override;
};