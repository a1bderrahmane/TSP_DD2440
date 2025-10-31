#include "Christofides.hpp"
#include "ITSPsolver.hpp"

using namespace std;

vector<int> Christofides::solve()
{
    return {};
}

Christofides::Christofides(const vector<Point> &points)
{
    int n = points.size();
    for (int i = 0; i < n; i++)
    {
        this->completeGraph.push_back({});
        for (int j = 0; j < n; j++)
        {
            this->completeGraph.back().push_back(computeDistance(points[i], points[j]));
        }
    }
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
        // Find vertex not yet in MST with smallest edge weight
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

        // Update edges to remaining vertices
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

    // Add edges based on parent relationships
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
