#include "GreedyTSPsolver.hpp"
#include "ITSPsolver.hpp"
#include "LKH.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
using namespace std;
using namespace chrono;

const double TIME_LIMIT = 1.80;    // seconds
const int CANDIDATE_SET_SIZE = 30; // number of nearest neighbors to consider

vector<int> LKH::solve(const vector<Point> &points)
{
    auto start = high_resolution_clock::now();

    int length = points.size();
    if (length == 0)
        return {};

    

    // Use GreedyTSPSolver to get an initial tour
    GreedyTSPSolver greedySolver;
    vector<int> tour = greedySolver.solve(points);
    double initialLength = computeTourLength(points, tour);

    if (length < 3)
        return tour;

    // Compute distance matrix
    vector<vector<int>> distance_matrix = computeDistanceMatrix(points, length);
    vector<vector<int>> candidate_sets = computeCandidateSets(distance_matrix, length);
    vector<int> position_array = computePositionArray(tour, length);

    double runtime = 0.0;
    while (runtime < TIME_LIMIT)
    {
        tour = lkStart(tour, distance_matrix, candidate_sets, position_array, length);
        position_array = computePositionArray(tour, length);

        runtime = duration<double>(high_resolution_clock::now() - start).count();
    }

    double finalLength = computeTourLength(points, tour);
    cout << "LKH improved tour from length " << initialLength << " to " << finalLength << endl;

    return tour;
}

vector<int> LKH::lkStart(vector<int> &tour, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                         vector<int> &position_array, const int length)
{
    vector<vector<int>> neighbors_matrix(length, vector<int>(2, -1));
    int city;
    for (int i = 0; i < length; i++)
    {
        city = tour[i];
        neighbors_matrix[city][0] = tour[(i - 1 + length) % length]; // backwards neighbor
        neighbors_matrix[city][1] = tour[(i + 1 + length) % length]; // forwards neighbor
    }
    MoveState best_move;
    best_move.neighbors = neighbors_matrix;

    for (int i = 0; i < length; i++)
    {
        int t1 = tour[i]; // pick a starting node

        // backward edge
        {
            MoveState backward;
            backward.neighbors = neighbors_matrix;
            int t2_pos = (i - 1 + length) % length;
            int t2 = tour[t2_pos];
            //backward.deletions.emplace_back(t1, t2);
            backward.remove(t1, t2);
            findNextY(tour, t1, t2, distanceMatrix, candidateSets, position_array, backward, best_move, length);
        }

        // forward edge
        {
            MoveState forward;
            forward.neighbors = neighbors_matrix;
            int t2_pos = (i + 1) % length;
            int t2 = tour[t2_pos];
            //forward.deletions.emplace_back(t1, t2);
            forward.remove(t1, t2);
            findNextY(tour, t1, t2, distanceMatrix, candidateSets, position_array, forward, best_move, length);
        }
    }

    if (best_move.gain > 0)
    {
        vector<int> best_tour;
        if (applyChanges(tour, best_tour, best_move, length))
        {
            tour = best_tour;
        }
    } // TODO: add random permutation if no gain is found to escape local minima

    return tour;
}

bool LKH::applyChanges(vector<int> &tour, vector<int> &best_tour, MoveState &best_move, int length)
{
    /* // getting neighbors in the current tour
    vector<vector<int>> neighbors_matrix(length, vector<int>(2, -1));
    int city;
    for (int i = 0; i < length; i++)
    {
        city = tour[i];
        neighbors_matrix[city][0] = tour[(i - 1 + length) % length]; // backwards neighbor
        neighbors_matrix[city][1] = tour[(i + 1 + length) % length]; // forwards neighbor
    }

    // applying the deletions
    for (int i = 0; i < best_move.deletions.size(); i++)
    {
        int a = best_move.deletions[i].first;
        int b = best_move.deletions[i].second;

        if (neighbors_matrix[a][0] == b)
        {
            neighbors_matrix[a][0] = -1;
        }
        else if (neighbors_matrix[a][1] == b)
        {
            neighbors_matrix[a][1] = -1;
        }
        else
        {
            return false;
        }

        if (neighbors_matrix[b][0] == a)
        {
            neighbors_matrix[b][0] = -1;
        }
        else if (neighbors_matrix[b][1] == a)
        {
            neighbors_matrix[b][1] = -1;
        }
        else
        {
            return false;
        }
    }

    // applying additions
    for (int i = 0; i < best_move.additions.size(); i++)
    {
        int a = best_move.additions[i].first;
        int b = best_move.additions[i].second;

        if (neighbors_matrix[a][0] == -1)
        {
            neighbors_matrix[a][0] = b;
        }
        else if (neighbors_matrix[a][1] == -1)
        {
            neighbors_matrix[a][1] = b;
        }
        else
        {
            return false;
        }

        if (neighbors_matrix[b][0] == -1)
        {
            neighbors_matrix[b][0] = a;
        }
        else if (neighbors_matrix[b][1] == -1)
        {
            neighbors_matrix[b][1] = a;
        }
        else
        {
            return false;
        }
    }
    */
    // rebuild the tour
    vector<vector<int>> neighbors_matrix = best_move.neighbors;

    cout << "Applying move with gain " << best_move.gain << endl;
    for (int i = 0; i < best_move.deletions.size(); i++)
    {
        cout << "  Deleting edge (" << best_move.deletions[i].first << ", " << best_move.deletions[i].second << ")" << endl;
    }
    for (int i = 0; i < best_move.additions.size(); i++)
    {
        cout << "  Adding edge (" << best_move.additions[i].first << ", " << best_move.additions[i].second << ")" << endl;
    }

    best_tour.assign(length, -1);
    int starting_node = tour[0];
    best_tour[0] = starting_node;
    int current_node = (neighbors_matrix[starting_node][0] != -1)
                           ? neighbors_matrix[starting_node][0]
                           : neighbors_matrix[starting_node][1];

    if (current_node == -1)
        return false; // broken adjacency

    best_tour[1] = current_node;

    int i = 2;
    int last_node = starting_node;
    int next_node;

    while (i < length)
    {
        if (neighbors_matrix[current_node][0] == last_node)
        {
            if (neighbors_matrix[current_node][1] != last_node)
            {
                next_node = neighbors_matrix[current_node][1];
            }
            else
            {
                return false;
            }
        }
        else if (neighbors_matrix[current_node][1] == last_node)
        {
            next_node = neighbors_matrix[current_node][0];
        }
        else
        {
            return false;
        }

        best_tour[i] = next_node;
        last_node = current_node;
        current_node = next_node;
        i += 1;
    }

    if (!isValidTour(best_tour, length))
    {
        return false;
    }
    cout << "Successfully applied move to obtain new tour." << endl;
    return true;
}

void LKH::findNextY(const vector<int> &tour, int t1, int current_vertex, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                    const vector<int> &position_array, MoveState &s, MoveState &best_move, int length)
{
    for (int i = 0; i < candidateSets[current_vertex].size(); i++)
    {
        int current_candidate = candidateSets[current_vertex][i];
        bool valid_candidate = validateYCandidate(tour, position_array, current_vertex, current_candidate, s, t1, length);

        if (!valid_candidate)
        {
            continue;
        }

        auto last_deleted_edge = s.deletions.back();
        double old_dist = distanceMatrix[last_deleted_edge.first][last_deleted_edge.second];
        double new_dist = distanceMatrix[current_vertex][candidateSets[current_vertex][i]];
        double local_gain = old_dist - new_dist;

        // TODO: we dont only want to allow positive gains (strict for y1, but later small negative ones are allowed to escape local minima)
        if (local_gain > 0) // this candidate gives a positive local gain
        {
            MoveState n = s;                                             // copy into a new candidate state
            n.gain += local_gain;                                        // update the gain
            //n.additions.emplace_back(current_vertex, current_candidate); // update the addition of the y edge
            n.add(current_vertex, current_candidate);
            findNextX(tour, t1, current_candidate, distanceMatrix, candidateSets, position_array, n, best_move, length);
        }
        else if (s.gain + local_gain > 0) // the tested candidate has a worse distance than the current link - all following candidates will be worse!
        {
            MoveState n = s;                                             // copy into a new candidate state
            n.gain += local_gain;                                        // update the gain
            // n.additions.emplace_back(current_vertex, current_candidate); // update the addition of the y edge
            n.add(current_vertex, current_candidate);
            findNextX(tour, t1, current_candidate, distanceMatrix, candidateSets, position_array, n, best_move, length);
        } else
        {
           break; // the loss is too high, all following candidates will be worse
        }
    }

    return; // we tried all candidates - if any found a better result, it will have been stored in best_move
}

void LKH::findNextX(const vector<int> &tour, int t1, int current_vertex, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                    const vector<int> &position_array, MoveState &s, MoveState &best_move, int length)
{
    int current_vertext_position = position_array[current_vertex];
    int forward_neighbor = tour[(current_vertext_position + 1) % length];
    int backward_neighbor = tour[(current_vertext_position - 1 + length) % length];

    int next_vertex = validateX(tour, t1, current_vertex, distanceMatrix, candidateSets, position_array, s, length, forward_neighbor, backward_neighbor, current_vertext_position);
    if (next_vertex == -1)
    { // no valid next edge to break found
        return;
    }

    // We found an edge (current_vertex, next_vertex) to delete

    // Calculate the total gain if we were to close now:
    double deleted_edge_distance = distanceMatrix[current_vertex][next_vertex];
    double closing_edge_distance = distanceMatrix[next_vertex][t1];
    double total_gain = s.gain + (deleted_edge_distance - closing_edge_distance);

    //s.deletions.emplace_back(current_vertex, next_vertex);
    s.remove(current_vertex, next_vertex);

    if (total_gain > best_move.gain) // in case we dont find any better solution we add full "instructions"
    {
        best_move = s;
        best_move.gain = total_gain;
        // best_move.additions.emplace_back(next_vertex, t1);
        best_move.add(next_vertex, t1);
    }

    if (s.additions.size() > 10) // limit the depth of the search to avoid too long runtimes
    {
        return;
    }

    // call findY function to get the next round of possible edge additions
    findNextY(tour, t1, next_vertex, distanceMatrix, candidateSets, position_array, s, best_move, length);
}

int LKH::validateX(const vector<int> &tour, int t1, int current_vertex, const vector<vector<int>> &distanceMatrix, const vector<vector<int>> &candidateSets,
                   const vector<int> &position_array, MoveState &s, int length, int forward_neighbor, int backward_neighbor, int t2i_minus_1_pos)
{
    int t2i_minus_1 = current_vertex;
    int t2i;
    int next_endpoint;

    if (!validateXCandidate(tour, position_array, current_vertex, forward_neighbor, s, length) || forward_neighbor == t1)
        {
            if (!validateXCandidate(tour, position_array, current_vertex, backward_neighbor, s, length) || backward_neighbor == t1)
            {
                return -1; // both forward and backward are not valid in x2
            }
            else
            { // backward edge is valid
                t2i = backward_neighbor;
                if (validateTourWithX(tour, position_array, t1, t2i, t2i_minus_1, s, length))
                {
                    return t2i; // a valid move has been found
                }
            }
        }
    
    if ((validateXCandidate(tour, position_array, current_vertex, forward_neighbor, s, length) || forward_neighbor != t1))
        { // forward edge is valid
            t2i = forward_neighbor;
            if (validateTourWithX(tour, position_array, t1, t2i, t2i_minus_1, s, length))
                {
                    return t2i; // a valid move has been found
                } else {
                    return -1;
                }
        } else {
            return -1;
        }


    
}; 


bool LKH::validateTourWithX(const vector<int> &tour, const vector<int> &position_array, int t1, int t2i,
                             int t2i_minus_1, const MoveState &s, int length) 
{
   
        int t2i_pos = position_array[t2i]; //FIXME check if we need to update the position array or if that works as is
        int t2i_minus_1_pos = position_array[t2i_minus_1];
        int stepDir;
        int i = t2i_minus_1_pos;

        if ((t2i_minus_1_pos + 1) % length == t2i_pos)
        {
            stepDir = -1; // t2i is forward -> walk backwards
        }
        else if ((t2i_minus_1_pos - 1 + length) % length == t2i_pos)
        {
            stepDir = +1; // t2i is backward -> walk forwards
        }
        else
        {
            return false; // t2i not neighbor of t2i_minus_1
        }

        int j = 0;
        int current_vertex = t2i_minus_1;
        int last_vertex = t2i;
        while (j < length) {
            j += 1;
            int current_vertex_neighbors_left = s.neighbors[current_vertex][0];
            int current_vertex_neighbors_right = s.neighbors[current_vertex][1];

            if (current_vertex_neighbors_left == last_vertex) {
                current_vertex = current_vertex_neighbors_right;
            } else if (current_vertex_neighbors_right == last_vertex) {
                current_vertex = current_vertex_neighbors_left;
            } else {
                return false; // broken adjacency
            }
            
            if (current_vertex == t1) {
                // we have an open edge here, no cycle
                return true;
            }

            if (current_vertex == t2i) {
                // we have a cycle before reaching t1
                return false;
            }
            
        } 
        return false; // we did not reach t1 within length steps 
};

// Validates a candidate for the closing Y edge:
bool LKH::validateYCandidate(const vector<int> &tour, const vector<int> &position_array, int current_vertex,
                             int current_candidate, const MoveState s, int t1, int length)
{
    // TODO - can be checked faster with hash map or sth like that - O(1) instead of O(n)
    auto edge_used = [&](const vector<pair<int,int>>& edges,
                     int a, int b) {
    for (const auto &p : edges) {
        if ((p.first == a && p.second == b) ||
            (p.first == b && p.second == a)) {
            return true;
        }
    }
    return false;
    };

    // we cannot close the loop too early - connect back to t1 only at the end
    if (current_candidate == t1)
    {
        return false;
    }

    if (edge_used(s.deletions, current_vertex, current_candidate))
        return false;

    if (edge_used(s.additions, current_vertex, current_candidate))
        return false;

    // edge is in the tour already
    int current_vertex_index = position_array[current_vertex];
    int left_neighbor = tour[(current_vertex_index - 1 + length) % length];
    int right_neighbor = tour[(current_vertex_index + 1) % length];
    if (current_candidate == left_neighbor || current_candidate == right_neighbor)
    {
        return false;
    }

    return true;
}

bool LKH::validateXCandidate(const vector<int> &tour, const vector<int> &position_array, int current_vertex,
                             int current_candidate, const MoveState &s, int length)
{
    // TODO: Does this need to be this strict or can we allow more edges?
    // check if the other end node of x has already been used, if yes, reject:
    // edge containing current_candidate has already been deleted
    for (const auto &p : s.deletions)
    {
        if (p.first == current_candidate || p.second == current_candidate)
        {
            return false;
        }
    }

    // edge containing current_candidate has already been added
    for (const auto &p : s.additions)
    {
        if (p.first == current_candidate || p.second == current_candidate)
        {
            return false;
        }
    }

    return true;
}

// Computes the full distance matrix between all points
vector<vector<int>> LKH::computeDistanceMatrix(const vector<Point> &points, int length)
{
    vector<vector<int>> dist(length, vector<int>(length, 0.0));

    for (int i = 0; i < length; i++)
    {
        for (int j = i + 1; j < length; j++)
        {
            int d = computeDistance(points[i], points[j]);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }
    return dist;
}

vector<int> LKH::computePositionArray(const vector<int> &tour, int length)
{
    vector<int> pos_array(length, 0);
    for (int i = 0; i < length; i++)
    {
        pos_array[tour[i]] = i;
    }
    return pos_array;
}

vector<vector<int>> LKH::computeCandidateSets(const vector<vector<int>> &distanceMatrix, int length)
{
    const int n = length;
    vector<vector<int>> candidate_sets(n);

    for (int i = 0; i < n; ++i)
    {
        // Build (distance, index) list for node i
        vector<pair<double, int>> distance_index_pairs;
        distance_index_pairs.reserve(n - 1);

        for (int j = 0; j < n; ++j)
        {
            if (i == j)
                continue;
            distance_index_pairs.emplace_back(distanceMatrix[i][j], j);
        }

        // Keep only the K closest neighbors
        if (CANDIDATE_SET_SIZE < length - 1)
        {
            nth_element(distance_index_pairs.begin(),
                        distance_index_pairs.begin() + CANDIDATE_SET_SIZE,
                        distance_index_pairs.end(),
                        [](const auto &a, const auto &b)
                        {
                            return a.first < b.first; // compare by distance
                        });
            distance_index_pairs.resize(CANDIDATE_SET_SIZE);
        }

        // Sort the selected neighbors by distance
        sort(distance_index_pairs.begin(), distance_index_pairs.end(),
             [](const auto &a, const auto &b)
             {
                 return a.first < b.first; // compare by distance
             });

        // Save the neighbor indices into candidate_sets[i]
        candidate_sets[i].reserve(CANDIDATE_SET_SIZE);
        for (const auto &pair : distance_index_pairs)
        {
            candidate_sets[i].push_back(pair.second);
        }
    }

    return candidate_sets;
}
// TODO
/*
1 - Create a pos_array lookup table: city 0 index in the tour should be stored at pos_array[0] CHECKED
2 - Compute a candidate set for each node (for example the k neighbors with shortest distance) CHECKED
3 - Start function which takes the tour and runs a for loop over all nodes to try out all pos_arraysible start nodes t1. It can also "delete"
(not actually but pick or mark) the edge between itself and the forward / backwards neighbor t2 and hands t2 to the next function (4).
->  Should also enforce time limit?? maybe a level below this function
->  Should consider which one of the starting nodes had the best run
->  Should create lists to keep track of removed and added edges
->  Should keep track of gains
4 - A function that finds a suitable y1 for a given node ti.

!! Change the name of the candidate handed down
!! The gain right now needs to be computed from t1-tn and t2-t3



5 - A function that finds the suitable xi+1 for a given yi (contacts C to check if needed).
C - A function that checks if a given move is valid or not.
*Functions 4 and 5 loop back and worth until one of them cant find a suitable node / edge anymore*
6 - A function that applies the found move to the tour.
7 - Update the start function to call function 3 until we cant find any improvement anymore or we run out of time.
*/


// Hash done tour to check if we alreadz have seen it - if yes break and restart with random tour
// Sort edges to start with the biggest edges in x1 !!
// First gain should be applied, then start at the same starting node that gave the best gain - still run through all the starting ndoes, just in a different order 