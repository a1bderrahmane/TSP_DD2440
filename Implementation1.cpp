#include "Implementation1.hpp"
#include "ITSPsolver.hpp"

using namespace std;

vector<int> Implementation1::solve( vector<Point> &points)
{
    int n = points.size();
    if (n == 0)
        return {};

    vector<int> tour(n);
    vector<bool> used(n, false);

    tour[0] = 0;
    used[0] = true;

    for (int i = 1; i < n; ++i)
    {
        int best = -1;
        int bestDist = -1;

        for (int j = 0; j < n; ++j)
        {
            if (!used[j])
            {
                int dist = computeDistance(points[tour[i - 1]], points[j]);
                if (best == -1 || dist < bestDist)
                {
                    best = j;
                    bestDist = dist;
                }
            }
        }

        tour[i] = best;
        used[best] = true;
    }

    return tour;
}
