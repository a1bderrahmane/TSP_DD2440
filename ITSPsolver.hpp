#ifndef TSP_SOLVER_INTERFACE_H
#define TSP_SOLVER_INTERFACE_H

#include <vector>
#include <cmath>
#include <utility>
#include <string>
using namespace std;

struct Point
{
    double x;
    double y;

    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
    bool operator==(const Point& other) const {
        return (this->x == other.x && this->y == other.y);
    }
};

// Computes Euclidean distance rounded to nearest integer
inline double computeDistance(const Point &a, const Point &b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

/**
 * Abstract base class for TSP solvers.
 * Implement this interface to create different solving strategies.
 */
class ITSPSolver
{
public:
    virtual ~ITSPSolver() = default;

    /**
     * Solves the TSP problem for the given points.
     * @param points Vector of 2D points to visit
     * @return Vector of indices representing the tour order (0-based)
     */
    virtual vector<int> solve(const vector<Point> &points) = 0;

    /**
     * Optional: Initialize solver with points (for preprocessing).
     * Default implementation does nothing.
     * @param points Vector of 2D points
     */
    virtual void initialize( vector<Point> &points)
    {
    }

protected:
    /**
     * Helper: Compute total tour length.
     * @param points Vector of all points
     * @param tour Vector of indices representing the tour
     * @return Total tour length
     */
     double computeTourLength(const vector<Point> &points,
                             const vector<int> &tour) const
    {
        if (tour.empty())
            return 0;

        double length = 0;
        for (size_t i = 0; i < tour.size() - 1; ++i)
        {
            length += computeDistance(points[tour[i]], points[tour[i + 1]]);
        }
        // Add distance from last point back to first
        length += computeDistance(points[tour.back()], points[tour[0]]);
        return length;
    }

    /**
     * Helper: Validate that a tour visits all points exactly once.
     * @param tour Vector of indices
     * @param n Number of points
     * @return true if valid, false otherwise
     */
    bool isValidTour(const vector<int> &tour, int n) const
    {
        if (tour.size() != n)
            return false;

        vector<bool> visited(n, false);
        for (int idx : tour)
        {
            if (idx < 0 || idx >= n || visited[idx])
                return false;
            visited[idx] = true;
        }
        return true;
    }
};

#endif
