#include "ITSPsolver.hpp"
#include <chrono>

struct MoveState
{
    double gain = 0.0;
    vector<pair<int, int>> additions;
    vector<pair<int, int>> deletions;
    vector<vector<int>> neighbors;

    void add(int a, int b) {
    additions.emplace_back(a, b);
    neighbors[a].push_back(b);
    neighbors[b].push_back(a);
    }

    void remove(int a, int b)
{
    // Record the deletion
    deletions.emplace_back(a, b);

    // ---- Remove b from neighbors[a] ----
    {
        bool found = false;
        for (size_t i = 0; i < neighbors[a].size(); i++)
        {
            if (neighbors[a][i] == b)
            {
                // Erase this entry
                neighbors[a].erase(neighbors[a].begin() + i);
                found = true;
                break;
            }
        }
    }

    // ---- Remove a from neighbors[b] ----
    {
        bool found = false;
        for (size_t i = 0; i < neighbors[b].size(); i++)
        {
            if (neighbors[b][i] == a)
            {
                neighbors[b].erase(neighbors[b].begin() + i);
                found = true;
                break;
            }
        }
    }
}

};

class LKH : public ITSPSolver
{
public:
    vector<int> solve(const vector<Point> &points) override;
    vector<int> lkStart(vector<int> &tour, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                        vector<int> &position_array, const int length);
    bool applyChanges(vector<int> &tour, vector<int> &new_tour, MoveState &best_move, int length);
    void findNextY(const vector<int> &tour, int starting_edge, int current_vertex, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                   const vector<int> &position_array, MoveState &s, MoveState &best_move, int length);
    void findNextX(const vector<int> &tour, int starting_edge, int current_vertex, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                   const vector<int> &position_array, MoveState &s, MoveState &best_move, int length);
    int validateX(const vector<int> &tour, int t1, int current_vertex, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                  const vector<int> &position_array, MoveState &s, int length, int forward_neighbor, int backward_neighbor, int current_vertex_position);
    bool validateYCandidate(const vector<int> &tour, const vector<int> &position_array, int current_vertex,
                            int current_candidate, const MoveState s, int t1, int length);
    bool validateXCandidate(const vector<int> &tour, const vector<int> &position_array, int current_vertex,
                            int current_candidate, const MoveState &s, int length);
    bool validateTourWithX(const vector<int> &tour, const vector<int> &position_array, int t1, int t2i,
                             int t2i_minus_1, const MoveState &s, int length);
    vector<vector<int>> computeDistanceMatrix(const vector<Point> &points, int length);
    vector<int> computePositionArray(const vector<int> &tour, int length);
    vector<vector<int>> computeCandidateSets(const vector<vector<int>> &distanceMatrix, int length);
};