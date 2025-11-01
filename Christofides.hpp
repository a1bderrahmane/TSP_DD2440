#include "ITSPsolver.hpp"

using namespace std;
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
struct MST
{
    vector<Edge> edges;
    void addEdge(Point u, Point v, double weight)
    {
        edges.push_back({u, v, weight});
    }
};


class Christofides : public ITSPSolver
{
public:
    Christofides(const vector<Point> &points);
    vector<int> solve(const vector<Point> &points) override;

private:
    MST mst;

    MST constructMST(const vector<Point> &points);
    vector<vector<double>> completeGraph;
    vector<int> getOddDegreeVertices(const vector<Point> &points, const MST &mst);
    vector<Edge> minWeightPerfectMatching(const vector<Point> &points, const vector<int> &oddVertices);
    vector<vector<pair<int, double>>> combineGraphs(const vector<Point> &points, const MST &mst, const vector<Edge> &matching);
    vector<int> findEulerianCircuit(vector<vector<pair<int, double>>> &graph);
    vector<int> makeHamiltonian(const vector<int> &eulerianCircuit);
};
