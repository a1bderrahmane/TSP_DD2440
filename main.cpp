#include "Christofides.hpp"
#include "Implementation2.hpp"
#include "Implementation3.hpp"
#include "Implementation4.hpp"
#include "GreedyTSPsolver.hpp"
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;


int main(int argc, char* argv[]) {
    // Check if implementation number is provided
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <implementation_number>" << endl;
        cerr << "  0 - GreedyTSPSolver (baseline)" << endl;
        cerr << "  1 - Christofides" << endl;
        cerr << "  2 - Implementation2" << endl;
        cerr << "  3 - Implementation3" << endl;
        cerr << "  4 - Implementation4" << endl;
        return 1;
    }
    
    // Parse implementation number
    int implNum = atoi(argv[1]);
    
    // Read number of points
    int n=2;
    
    
    // Read all points
    vector<Point> points(n);
    points[0].x=1;
    points[0].y=2;
    points[1].x=3;
    points[1].y=4;
    
    // Solve based on implementation number
    vector<int> tour;
    
    switch (implNum) {
        case 0: {
            GreedyTSPSolver solver;
            tour = solver.solve(points);
            break;
        }
        case 1: {
            Christofides solver= Christofides(points);
            tour = solver.solve(points);
            break;
        }
        case 2: {
            // Implementation2 solver;
            // tour = solver.solve(points);
            break;
        }
        case 3: {
            // Implementation3 solver;
            // tour = solver.solve(points);
            break;
        }
        case 4: {
            // Implementation4 solver;
            // tour = solver.solve(points);
            break;
        }
        default:
            cerr << "Error: Invalid implementation number " << implNum << endl;
            cerr << "Valid options: 0 (Greedy), 1, 2, 3, 4" << endl;
            return 1;
    }
    
    // Output the tour (one index per line)
    for (int idx : tour) {
        cout << idx << "\n";
    }
    
    return 0;
}