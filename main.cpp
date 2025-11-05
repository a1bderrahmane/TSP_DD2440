#include "LKH.hpp"
#include <iostream>
#include <vector>
using namespace std;

int main()
{
    int n;
    if (!(cin >> n))
    {
        return 0; // no input -> just exit
    }

    vector<Point> points(n);
    for (int i = 0; i < n; ++i)
    {
        cin >> points[i].x >> points[i].y;
    }

    LKH solver;
    vector<int> tour = solver.solve(points);

    for (int idx : tour)
    {
        cout << idx << "\n";
    }

    return 0;
}
