#include "ITSPsolver.hpp"

using namespace std;

class Implementation1 : public ITSPSolver
{
public:
    vector<int> solve(vector<Point> &points) override;
};
