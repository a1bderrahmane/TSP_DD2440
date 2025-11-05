#include "GreedyTSPsolver.hpp"
#include "ITSPsolver.hpp"
#include "LKH.hpp"
#include <algorithm>

vector<int> LKH::solve(const vector<Point> &points)
{

    int n = points.size();
    if (n == 0)
        return {};

    // Use GreedyTSPSolver to get an initial tour
    GreedyTSPSolver greedySolver;
    vector<int> tour = greedySolver.solve(points);

    bool improved = true;
    int counter = 1; // counter to prevent running out of time

    while (improved)
    {
        improved = false;
        // check all pairs of edges for possible 2-opt swaps
        for (int i = 1; i < n - 1; ++i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                int a = tour[i - 1];
                int b = tour[i];
                int c = tour[j];
                int d = tour[(j + 1) % n];

                int currentDistance = computeDistance(points[a], points[b]) +
                                      computeDistance(points[c], points[d]);

                int newDistance = computeDistance(points[a], points[c]) +
                                  computeDistance(points[b], points[d]);

                // If the new distance is shorter, perform the 2-opt swap by reversing the segment between i and j
                if (newDistance < currentDistance)
                {
                    std::reverse(tour.begin() + i, tour.begin() + j + 1);
                    improved = true;
                }
            }
        }
        if (counter > 200)
            break;
        counter++;
    }

    return tour;
}