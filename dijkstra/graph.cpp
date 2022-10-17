/* Graph Class for Dijkstra's Algorithm
 * 
 * Week 3 assignment - C++ for C programmers
 * 
 */


#include <random>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include "graph.h"


using namespace dijkstra;


MinHeap::MinHeap(): queue(compare_second), val_map()
{}


void MinHeap::set(int index, double value)
{
    // Note we don't go through the priority queue to remove the outdated value.
    val_map[index] = value;
    queue.push(queue_val(index, value));
}


double MinHeap::get(int index) const
{
    // Note unordered_map operator[] is not const because it will add an element if it is not found.
    return val_map.at(index);
}


double MinHeap::operator[](int index) const
{
    // Note this doesn't allow setting.  The set version of this return a reference to the data.  The value added to
    // the priority_queue shouldn't be changed after pushing in.  To avoid this, just don't allow setting and make
    // const.
    return get(index);
}


double MinHeap::contains(int index) const
{
    return val_map.find(index) != val_map.end();
}


void MinHeap::erase(int index)
{
    val_map.erase(index);
}


queue_val MinHeap::get_min() {
    queue_val out;
    while (!queue.empty()) {
        out = queue.top();
        queue.pop();
        // Instead of removing oudated values when setting, just skip outdated weights that don't match values in the
        // val_map.  This avoids rebalancing the priority_queue.
        if (val_map.find(out.first) != val_map.end() && val_map[out.first] == out.second) {
            val_map.erase(out.first);
            return out;
        }
    }
    // If we're here the queue doesn't have any useful elements.
    out = {-1, 0.0};
    return out;
}


int MinHeap::size() const
{
    return val_map.size();
}


Vertex::Vertex(int id): id(id), conn_list(0, nullptr)
{}


int Vertex::get_id() const
{
    return id;
}


void Vertex::add_edge(std::shared_ptr<Edge> edge)
{
    conn_list.push_back(edge);
}


const std::vector<std::pair<int, double>> Vertex::get_edges() const
{
    // We probably could precompute this for the version of the graphs we have.  However, I'm going to assume users
    // can add new nodes to the graph, so precomputing isn't the best idea.
    std::vector<std::pair<int, double>> out_list;
    for (auto edge : conn_list) {
        int other_i = edge->get_other(*this).get_id();
        double weight = edge->get_weight();
        out_list.push_back(std::pair<int, double>(other_i, weight));
    }
    return out_list;
}


Edge::Edge(std::shared_ptr<Vertex> vertex_a, std::shared_ptr<Vertex> vertex_b, double weight):
    vertex_a(vertex_a), vertex_b(vertex_b), weight(weight)
{}


const Vertex &Edge::get_a() const
{
    return *vertex_a;
}


const Vertex &Edge::get_b() const
{
    return *vertex_b;
}


const Vertex &Edge::get_other(const Vertex &thing) const
{
    if(thing.get_id() == vertex_a->get_id()) {
        return get_b();
    }
    else if (thing.get_id() == vertex_b->get_id()) {
        return get_a();
    }
    else {
        throw std::invalid_argument("Input vertex not in edge");
    }
}


double Edge::get_weight() const
{
    return weight;
}


Path::Path(double weight): weight(weight), path()
{}


double Path::get_weight() const
{
    return weight;
}


void Path::push_front(int id)
{
    path.push_front(id);
}


void Path::push_back(int id)
{
    path.push_back(id);
}


const std::list<int> Path::get_path() const
{
    // all of the "const" on here should ensure this is immutable?
    return path;
}


Graph::Graph(int n_nodes, double density, double dist_min, double dist_max):
    n_nodes(n_nodes), density(density), dist_min(dist_min), dist_max(dist_max), vertex_list()
{
    // Initialize random number generator for the graph
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> edge_dist(0.0, 1.0);
    std::uniform_real_distribution<double> distance_dist(dist_min, dist_max);

    // Generate Node List
    for (int node_i=0; node_i<n_nodes; ++node_i) {
        vertex_list.push_back(std::make_shared<Vertex>(node_i));
    }

    // Generate Edges
    for (int cur_i=0; cur_i<n_nodes-1; ++cur_i) {
        for (int other_i=cur_i+1; other_i<n_nodes; ++other_i) {
            if(edge_dist(generator) > density) {
                continue;
            }
            double edge_weight = distance_dist(generator);
            auto cur_edge = std::make_shared<Edge>(
                vertex_list[cur_i],
                vertex_list[other_i],
                edge_weight);
            vertex_list[cur_i]->add_edge(cur_edge);
            vertex_list[other_i]->add_edge(cur_edge);
        }
    }
}


const Path Graph::find_path(int start, int end) const
{
    int cur_i = start;
    double cur_weight = 0.0;
    std::unordered_set<int> closed_set = {};  // only need to remember which nodes have been visited
    MinHeap open_weight;  // Min heap to make finding minimum value fast
    std::unordered_map<int, int> prev_dict = {};  // Points to previous node

    open_weight.set(start, 0.0);

    // Limit the number of nodes we check to the number of nodes, in case I do something silly
    for (int iter_i=0; iter_i < vertex_list.size(); ++iter_i) {
        // We're done if we're at the end
        if (cur_i == end) {
            break;
        }
        auto out_edges = vertex_list[cur_i]->get_edges();
        // Check each edge leaving the vertex
        for (auto edge: out_edges) {
            // Don't bother checking if we've already closed the other side of the edge
            if(closed_set.find(edge.first) != closed_set.end()) {
                continue;
            }
            double new_weight = cur_weight + edge.second;
            if (open_weight.contains(edge.first)) {
                // vertex has already been visited, update if this has a lower cost than previous visit.
                if(new_weight < open_weight[edge.first]) {
                    open_weight.set(edge.first, new_weight);
                    prev_dict[edge.first] = cur_i;
                }
            }
            else {
                // Never visited this vertex, so update
                open_weight.set(edge.first, new_weight);
                prev_dict[edge.first] = cur_i;
            }
        }
        // Done with current vertex, add to closed set and remove from possible nodes.
        closed_set.insert(cur_i);
        open_weight.erase(cur_i);
        // Find next vertex to work from (lowest weight thus far)
        std::tie(cur_i, cur_weight) = open_weight.get_min();
        // No more nodes to check, so end must not be reachable from start.
        if (cur_i == -1) {
            // Not sure how to do exceptions yet, so I'm just going to return an empty path
            std::cout << "No Path found!" << std::endl;
            return Path(-1);
        }
    }
    // If we're here, we found a path
    Path out_path(cur_weight);
    for (int iter_i=0; iter_i < vertex_list.size(); ++iter_i) {
        out_path.push_front(cur_i);
        if (cur_i == start) {
            break;
        }
        cur_i = prev_dict[cur_i];
    }
    return out_path;
}


/**
 * Shortest Path Test
*/
class ShortestPathTest {
    public:
    int n_graphs;
    int n_paths;
    int n_nodes;
    double density;
    double dist_min;
    double dist_max;
    bool random_paths;
    bool debug_print;

    ShortestPathTest(int n_graphs, int n_paths, int n_nodes, double density, double dist_min, double dist_max,
                     bool random_paths, bool debug_print):
            n_graphs(n_graphs), n_paths(n_paths), n_nodes(n_nodes), density(density), dist_min(dist_min),
            dist_max(dist_max), random_paths(random_paths), debug_print(debug_print)
    {}

    friend std::ostream &operator<<(std::ostream &os, ShortestPathTest const &setup)
    {
        os << "GraphTestSetup(n_graphs=" << setup.n_graphs
            << ", n_paths=" << setup.n_paths
            << ", n_nodes=" << setup.n_nodes
            << ", density=" << setup.density
            << ", dist_min=" << setup.dist_min << ", dist_max=" << setup.dist_max
            << ", random_paths=" << setup.random_paths
            << ", debug_print=" << setup.debug_print << ")";
        return os;
    };

    double test_shortest_path()
    {
        std::cout << "Running " << *this << "\n";
        // Initialize path setup.  Pretty sure this grabs system clock or something to initialize, so it should be fine
        // setting up the same way as above.
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_int_distribution<int> path_dist(0, n_nodes-1);  // Dist is [0, n-1] inclusive

        int n_paths = n_nodes * (n_nodes - 1) / 2;
        std::vector<std::pair<int, int>> path_pair;
        if (random_paths) {
            n_paths = n_paths;
        }
        else {
            // Number of paths between nodes is the number of nodes squared (all pairs of nodes), without paths from
            // nodes back to itself (n_nodes) and symmetric paths (i to j == j to i).  This is the combinatorics thing
            // of N-choose-2.
            n_paths = n_nodes * (n_nodes - 1) / 2;
            for (int path_a = 0; path_a < n_nodes-1; ++path_a) {
                for (int path_b = path_a + 1; path_b < n_nodes; ++path_b) {
                    path_pair.push_back(std::pair<int, int>(path_a, path_b));
                }
            }
        }

        // Generate paths
        int n_trials = n_graphs * n_paths;
        std::vector<double> mean_weight_list;  // mean weight for each graph
        for (int graph_i = 0; graph_i < n_graphs; ++graph_i) {
            int n_paths_found = 0;
            double total_weight = 0.0;
            // Initialize graph.  Should be freed up after going out of scope (each iteration).
            Graph graph(n_nodes, density, dist_min, dist_max);
            if (debug_print) {
                std::cout << "Graph " << graph_i << "\n" << graph;
            }
            // Iterate over paths, sampling randomly or doing every possible.
            for (int path_i = 0; path_i < n_paths; ++path_i) {
                int start, end;
                if (random_paths) {
                    start = path_dist(generator);
                    end = path_dist(generator);
                    while (start == end) {
                        end = path_dist(generator);
                    }
                }
                else {
                    std::tie(start, end) = path_pair[path_i];
                }
                if (debug_print) {
                    std::cout << "Find Path from " << start << " to " << end << "\n";
                }
                Path out_path = graph.find_path(start, end);
                double weight = out_path.get_weight();
                if (weight > 0) {
                    if (debug_print) {
                        std::cout << "Path From " << start << " to " << end << " Found\n" << out_path;
                    }
                    ++n_paths_found;
                    total_weight += weight;
                }
                else if (debug_print) {
                    std::cout << "No Path from " << start << " to " << end << " Found\n";
                }
            }
            mean_weight_list.push_back(total_weight / n_paths_found);
        }
        // Compute mean weight of all graphs and standard deviation of weight means (Central Limit Theorem, etc etc)
        double weight_mean = std::accumulate(mean_weight_list.begin(), mean_weight_list.end(), 0.0);
        weight_mean /= mean_weight_list.size();
        double weight_var = 0.0;
        for (auto mw: mean_weight_list) {
            weight_var += std::pow(mw - weight_mean, 2);
        }
        weight_var /= mean_weight_list.size() - 1;
        double weight_std = std::sqrt(weight_var);
        std::cout << " - Mean Weight = " << weight_mean << "\n - Std. Dev.   = " << weight_std << std::endl;
        return weight_mean;
    }
};


int main() {
    // Main setup parameters
    ShortestPathTest setup(
        100,  // n_graphs
        10,  // n_paths
        50,  // n_nodes
        0.4,  // density
        1.0,  // dist_min
        10.0,  // dist_max
        false,  // random_paths
        false // debug_print
    );
    bool run_output = true;
    if (run_output) {
        // Run for homework results (with different densities)
        setup.density = 0.4;
        setup.test_shortest_path();
        setup.density = 0.2;
        setup.test_shortest_path();
    }
    else {
        // Run debugging setup
        setup.test_shortest_path();
    }
    return 0;
}


