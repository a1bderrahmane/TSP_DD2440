#include "Christofides.hpp"
#include <limits>
#include <cmath>
#include <unordered_map>
#include <set>
#include <stack>
#include <algorithm>

using namespace std;

Christofides::Christofides(const vector<Point> &points)
{
    int n = points.size();
    this->completeGraph.resize(n, vector<double>(n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            completeGraph[i][j] = computeDistance(points[i], points[j]);
        }
    }
}

vector<int> Christofides::solve(const vector<Point> &points)
{
    // Construct MST
    mst = constructMST(points);

    // find nodes with odd degree
    vector<int> oddVertices = getOddDegreeVertices(points, mst);

    // Find minimum weight perfect matching among odd degree vertices
    vector<Edge> matching = minWeightPerfectMatching(points, oddVertices);

    // Combine MST + matching to form Eulerian multigraph
    vector<vector<pair<int, double>>> eulerianGraph = combineGraphs(points, mst, matching);

    // Find Eulerian circuit
    vector<int> eulerianCircuit = findEulerianCircuit(eulerianGraph);

    // Shortcut repeated nodes to form TSP tour
    vector<int> tspTour = makeHamiltonian(eulerianCircuit);

    return tspTour;
}

MST Christofides::constructMST(const vector<Point> &points)
{
    int n = completeGraph.size();
    vector<double> minEdge(n, numeric_limits<double>::infinity());
    vector<int> parent(n, -1);
    vector<bool> inMST(n, false);

    minEdge[0] = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double minWeight = numeric_limits<double>::infinity();
        int u = -1;

        for (int v = 0; v < n; ++v)
        {
            if (!inMST[v] && minEdge[v] < minWeight)
            {
                minWeight = minEdge[v];
                u = v;
            }
        }

        if (u == -1)
            break;

        inMST[u] = true;

        for (int v = 0; v < n; ++v)
        {
            double w = completeGraph[u][v];
            if (!inMST[v] && w < minEdge[v])
            {
                minEdge[v] = w;
                parent[v] = u;
            }
        }
    }

    MST mst;
    for (int v = 1; v < n; ++v)
    {
        if (parent[v] != -1)
        {
            Point uPoint = points[parent[v]];
            Point vPoint = points[v];
            double w = completeGraph[parent[v]][v];
            mst.addEdge(uPoint, vPoint, w);
        }
    }
    return mst;
}

vector<int> Christofides::getOddDegreeVertices(const vector<Point> &points, const MST &mst)
{
    unordered_map<int, int> degree;
    for (int i = 0; i < points.size(); ++i)
        degree[i] = 0;

    for (const auto &e : mst.edges)
    {
        int i1 = find(points.begin(), points.end(), e.node1) - points.begin();
        int i2 = find(points.begin(), points.end(), e.node2) - points.begin();
        degree[i1]++;
        degree[i2]++;
    }

    vector<int> oddVertices;
    for (auto &[idx, deg] : degree)
    {
        if (deg % 2 == 1)
            oddVertices.push_back(idx);
    }

    return oddVertices;
}

vector<Edge> Christofides::minWeightPerfectMatching(const vector<Point> &points, const vector<int> &oddVertices)
{
    vector<Edge> matching;
    set<int> unmatched(oddVertices.begin(), oddVertices.end());

    while (!unmatched.empty())
    {
        int u = *unmatched.begin();
        unmatched.erase(u);

        double bestWeight = numeric_limits<double>::infinity();
        int bestV = -1;

        for (int v : unmatched)
        {
            double w = completeGraph[u][v];
            if (w < bestWeight)
            {
                bestWeight = w;
                bestV = v;
            }
        }

        if (bestV != -1)
        {
            unmatched.erase(bestV);
            matching.push_back({points[u], points[bestV], bestWeight});
        }
    }
    return matching;
}

vector<vector<pair<int, double>>> Christofides::combineGraphs(const vector<Point> &points, const MST &mst, const vector<Edge> &matching)
{
    int n = points.size();
    vector<vector<pair<int, double>>> graph(n);

    auto addEdge = [&](int u, int v, double w)
    {
        graph[u].push_back({v, w});
        graph[v].push_back({u, w});
    };

    for (auto &e : mst.edges)
    {
        int u = find(points.begin(), points.end(), e.node1) - points.begin();
        int v = find(points.begin(), points.end(), e.node2) - points.begin();
        addEdge(u, v, e.weight);
    }

    for (auto &e : matching)
    {
        int u = find(points.begin(), points.end(), e.node1) - points.begin();
        int v = find(points.begin(), points.end(), e.node2) - points.begin();
        addEdge(u, v, e.weight);
    }

    return graph;
}

vector<int> Christofides::findEulerianCircuit(vector<vector<pair<int, double>>> &graph)
{
    vector<int> circuit;
    stack<int> currPath;
    currPath.push(0);
    vector<int> currentDegree(graph.size());
    for (int i = 0; i < graph.size(); ++i)
        currentDegree[i] = graph[i].size();

    while (!currPath.empty())
    {
        int u = currPath.top();
        if (currentDegree[u])
        {
            auto [v, w] = graph[u].back();
            graph[u].pop_back();
            currentDegree[u]--;
            currPath.push(v);
        }
        else
        {
            circuit.push_back(u);
            currPath.pop();
        }
    }

    reverse(circuit.begin(), circuit.end());
    return circuit;
}

vector<int> Christofides::makeHamiltonian(const vector<int> &eulerianCircuit)
{
    vector<int> path;
    set<int> visited;

    for (int v : eulerianCircuit)
    {
        if (visited.insert(v).second)
            path.push_back(v);
    }
    path.push_back(path[0]);
    return path;
}
