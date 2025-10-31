#include "ITSPsolver.hpp"

using namespace std;
struct MST
{
    vector<Edge> edges;
};

struct Edge
{
    Point node1;
    Point node2;
    bool operator==(const Edge &other) const
    {
        return (node1 == other.node1 && node2 == other.node2) ||
               (node1 == other.node2 && node2 == other.node1);
    }
};

class Christofides : public ITSPSolver
{
public:
    Christofides(const vector<Point> &points);
    vector<int> solve() override;

private:
    MST constructMST();
    vector<vector<double>> completeGraph;
    MST mst;
};
