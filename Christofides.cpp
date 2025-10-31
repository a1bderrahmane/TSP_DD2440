#include "Christofides.hpp"
#include "ITSPsolver.hpp"

using namespace std;

vector<int> Christofides::solve()
{
    return {};
}

Christofides::Christofides(const vector<Point> & points)
{
    int n =points.size();
    for(int i=0;i<n;i++)
    {
        this->completeGraph.push_back({});
        for(int j=0;j<n;j++)
        {
            this->completeGraph.back().push_back(computeDistance(points[i],points[j]));
        }
    }
}