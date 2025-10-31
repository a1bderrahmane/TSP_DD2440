#include "ITSPsolver.hpp"

using namespace std;

class Christofides : public ITSPSolver
{
public:
    vector<int> solve(vector<Point> &points) override;
};
