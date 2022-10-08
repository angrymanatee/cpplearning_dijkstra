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
    std::unordered_set<int> closed_set = {};
    std::unordered_map<int, double> open_weight = {};
    std::unordered_map<int, int> prev_dict = {};

    open_weight[start] = 0.0;

    // Limit the number of nodes we check to the number of nodes, in case I do something
    // silly
    for (int iter_i=0; iter_i < vertex_list.size(); ++iter_i) {
        // We're done if we're at the end
        if (cur_i == end) {
            break;
        }
        // Get the current weight and the edges out of this vertex
        double cur_weight = open_weight[cur_i];
        auto out_edges = vertex_list[cur_i]->get_edges();
        // Check each edge leaving the vertex
        for (auto edge: out_edges) {
            // Don't bother checking if we've already closed the other side of the edge
            if(closed_set.find(edge.first) != closed_set.end()) {
                continue;
            }
            double new_weight = cur_weight + edge.second;
            if (open_weight.find(edge.first) != open_weight.end()) {
                // vertex has already been visited, update if this has a lower cost than previous visit.
                if(new_weight < open_weight[edge.first]) {
                    open_weight[edge.first] = new_weight;
                    prev_dict[edge.first] = cur_i;
                }
            }
            else {
                // Never visited this vertex, so update
                open_weight[edge.first] = new_weight;
                prev_dict[edge.first] = cur_i;
            }
        }
        // Done with current vertex, add to closed set and remove from possible nodes.
        closed_set.insert(cur_i);
        open_weight.erase(cur_i);
        // Find next vertex to work from (lowest weight thus far)
        std::pair<int, double> cur_min_pair = {-1, INFINITY};
        for (auto test_weight_pair: open_weight) {
            // Checking closed should be unnecessary?
            // if (closed_set.find(test_weight_pair.first) != closed_set.end()) {
            //     continue;
            // }
            if (test_weight_pair.second < cur_min_pair.second) {
                cur_min_pair = test_weight_pair;
            }
        }
        // No more nodes to check, so end must not be reachable from start.
        cur_i = cur_min_pair.first;
        if (cur_i == -1) {
            // Not sure how to do exceptions yet, so I'm just going to return an empty path
            std::cout << "No Path found!" << std::endl;
            return Path(-1);
        }
    }
    // If we're here, we found a path
    Path out_path(open_weight[end]);
    cur_i = end;
    for (int iter_i=0; iter_i < vertex_list.size(); ++iter_i) {
        out_path.push_front(cur_i);
        if (cur_i == start) {
            break;
        }
        cur_i = prev_dict[cur_i];
    }
    return out_path;
}


struct TestSetup {
    int n_graphs;
    int n_paths;
    int n_nodes;
    double density;
    double dist_min;
    double dist_max;
    bool random_paths;
    bool debug_print;
    friend std::ostream &operator<<(std::ostream &os, TestSetup const &setup)
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
};


double run_test(const TestSetup &setup)
{
    std::cout << "Running " << setup << "\n";
    // Initialize path setup.  Pretty sure this grabs system clock or something to initialize, so it should be fine
    // setting up the same way as above.
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<int> path_dist(0, setup.n_nodes-1);  // Dist is [0, n-1] inclusive

    int n_paths = setup.n_nodes * (setup.n_nodes - 1) / 2;
    std::vector<std::pair<int, int>> path_pair;
    if (setup.random_paths) {
        n_paths = setup.n_paths;
    }
    else {
        // Number of paths between nodes is the number of nodes squared (all pairs of nodes), without paths from nodes
        // back to itself (n_nodes) and symmetric paths (i to j == j to i).  This is the combinatorics thing of
        // N-choose-2.
        n_paths = setup.n_nodes * (setup.n_nodes - 1) / 2;
        for (int path_a = 0; path_a < setup.n_nodes-1; ++path_a) {
            for (int path_b = path_a + 1; path_b < setup.n_nodes; ++path_b) {
                path_pair.push_back(std::pair<int, int>(path_a, path_b));
            }
        }
    }

    int n_trials = setup.n_graphs * n_paths;
    int n_paths_found = 0;
    double total_weight = 0.0;
    std::vector<double> weight_list(n_trials);
    for (int graph_i = 0; graph_i < setup.n_graphs; ++graph_i) {
        Graph graph(setup.n_nodes, setup.density, setup.dist_min, setup.dist_max);
        if (setup.debug_print) {
            std::cout << "Graph " << graph_i << "\n" << graph;
        }
        // Dump graph to screen or something for debugging?
        for (int path_i = 0; path_i < n_paths; ++path_i) {
            int start, end;
            if (setup.random_paths) {
                start = path_dist(generator);
                end = path_dist(generator);
                while (start == end) {
                    end = path_dist(generator);
                }
            }
            else {
                std::tie(start, end) = path_pair[path_i];
            }
            if (setup.debug_print) {
                std::cout << "Find Path from " << start << " to " << end << "\n";
            }
            Path out_path = graph.find_path(start, end);
            double weight = out_path.get_weight();
            if (weight > 0) {
                if (setup.debug_print) {
                    std::cout << "Path From " << start << " to " << end << " Found\n" << out_path;
                }
                ++n_paths_found;
                total_weight += weight;
            }
            else if (setup.debug_print) {
                std::cout << "No Path from " << start << " to " << end << " Found\n";
            }
        }
    }
    double weight_mean = total_weight / n_paths_found;
    std::cout << "Mean Weight = " << weight_mean << std::endl;
    return weight_mean;
}


int main() {
    TestSetup setup = {
        .n_graphs = 100,
        .n_paths = 10,
        .n_nodes = 50,
        .density = 0.4,
        .dist_min = 1.0,
        .dist_max = 10.0,
        .random_paths = false,
        .debug_print = false,
    };
    bool run_output = true;
    if (run_output) {
        setup.density = 0.4;
        std::cout << "Run Test, density = " << setup.density << std::endl;
        run_test(setup);

        setup.density = 0.2;
        std::cout << "Run Test, density = " << setup.density << std::endl;
        run_test(setup);
    }
    else {
        std::cout << "Run Test, density = " << setup.density << std::endl;
        run_test(setup);
    }
    return 0;
}
