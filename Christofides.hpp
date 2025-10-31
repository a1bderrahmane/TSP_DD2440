#include "ITSPsolver.hpp"

using namespace std;
struct MST
{
    vector<Edge> edges;
    void addEdge(Point u, Point v, double weight) {
        edges.push_back({u, v, weight});
    }
};

struct Edge
{
    Point node1;
    Point node2;
    double weight;
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
    MST constructMST(const vector<Point>&points);
    vector<vector<double>> completeGraph;
    MST mst;
};
