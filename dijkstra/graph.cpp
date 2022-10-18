/* Graph Class for Dijkstra's Algorithm
 * 
 * Week 3 assignment - C++ for C programmers
 * 
 */


#include <random>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <fstream>
#include <regex>
#include "graph.h"


using namespace dijkstra;


template <class T>
MinHeap<T>::MinHeap(): queue(compare_second<T>), val_map()
{}


template <class T>
void MinHeap<T>::set(T index, double value)
{
    // Note we don't go through the priority queue to remove the outdated value.
    val_map[index] = value;
    queue.push(MinHeap<T>::queue_val(index, value));
}


template <class T>
double MinHeap<T>::get(T index) const
{
    // Note unordered_map operator[] is not const because it will add an element if it is not found.
    return val_map.at(index);
}


template <class T>
double MinHeap<T>::operator[](T index) const
{
    // Note this doesn't allow setting.  The set version of this return a reference to the data.  The value added to
    // the priority_queue shouldn't be changed after pushing in.  To avoid this, just don't allow setting and make
    // const.
    return get(index);
}


template <class T>
double MinHeap<T>::contains(T index) const
{
    return val_map.find(index) != val_map.end();
}


template <class T>
void MinHeap<T>::erase(T index)
{
    val_map.erase(index);
}


template <class T>
typename MinHeap<T>::queue_val MinHeap<T>::get_min()
{
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


template <class T>
int MinHeap<T>::size() const
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


bool Vertex::is_connected_to(int id)
{
    for (auto edge : conn_list) {
        int other_i = edge->get_other(*this).get_id();
        if (other_i == id) {
            return true;
        }
    }
    return false;
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


Tree::Tree(): head(nullptr), vertex_map(), edge_list()
{}


Tree::Tree(int head_id): Tree()
{
    set_head(head_id);
}


void Tree::set_head(int head_id)
{
    if (head != nullptr) {
        throw std::invalid_argument("Head node already set");
    }
    head = std::make_shared<Vertex>(head_id);
    vertex_map[head_id] = head;
}


const Vertex &Tree::get_head() const
{
    return *head;
}


void Tree::add_child(int parent_id, int child_id, double weight)
{
    // Make sure child is not already a vertex
    if (vertex_map.find(child_id) != vertex_map.end()) {
        throw std::invalid_argument("Child already added");
    }
    // Create edge to child
    vertex_map[child_id] = std::make_shared<Vertex>(child_id);
    auto cur_edge = std::make_shared<Edge>(vertex_map[parent_id], vertex_map[child_id], weight);
    vertex_map[parent_id]->add_edge(cur_edge);
    vertex_map[child_id]->add_edge(cur_edge);
    edge_list.push_back(std::move(cur_edge));
}


double Tree::get_total_weight() const
{
    double total;
    for (auto edge_ptr : edge_list) {
        total += edge_ptr->get_weight();
    }
    return total;
}


void Tree::dump_children(std::ostream &os, int node_i, double weight, int n_indent) const
{
    for (int indent_i = 0; indent_i < n_indent; ++indent_i) {
        os << "|   ";
    }
    os << "|-(" << weight << ")-> " << node_i << "\n";
    for (auto edge_info : vertex_map.at(node_i)->get_edges()) {
        dump_children(os, edge_info.first, edge_info.second, n_indent+1);
    }
}


Graph::Graph(int n_nodes): n_nodes(n_nodes), vertex_list(), edge_list()
{
    create_node_list(n_nodes);
}


Graph::Graph(int n_nodes, double density, double dist_min, double dist_max): Graph(n_nodes)
{
    // Initialize random number generator for the graph
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> edge_dist(0.0, 1.0);
    std::uniform_real_distribution<double> distance_dist(dist_min, dist_max);

    // Generate Edges
    for (int cur_i=0; cur_i<n_nodes-1; ++cur_i) {
        for (int other_i=cur_i+1; other_i<n_nodes; ++other_i) {
            if(edge_dist(generator) > density) {
                continue;
            }
            double edge_weight = distance_dist(generator);
            add_edge(cur_i, other_i, edge_weight);
        }
    }
}


Graph::Graph(std::string graph_fname): n_nodes(), vertex_list(), edge_list()
{
    int max_char = 32;
    std::string line;
    std::smatch match;
    std::regex edge_re("(\\d+) (\\d+) (\\d+)\\s*");
    std::ifstream graph_file(graph_fname);

    // first line is the number of nodes
    std::getline(graph_file, line);
    n_nodes = std::stoi(line);
    create_node_list(n_nodes);

    // Rest of file is edges
    while (std::getline(graph_file, line)) {
        if (!std::regex_match(line, match, edge_re)) {
            continue;
        }
        // Note match[0] is full input string (regex thing)
        int node_a = std::stoi(match[1]);
        int node_b = std::stoi(match[2]);
        double weight = std::stod(match[3]);
        // Looks like the example file has redundant edges (a -> b and b -> a exist).  Since we know this is
        // undirected and symmetric, just ignore repeated edges.
        if (vertex_list[node_a]->is_connected_to(node_b)) {
            continue;
        }
        add_edge(node_a, node_b, weight);
    }
}


void Graph::create_node_list(int n_nodes)
{
    for (int node_i=0; node_i<n_nodes; ++node_i) {
        vertex_list.push_back(std::make_shared<Vertex>(node_i));
    }
}


void Graph::add_edge(int node_a, int node_b, double weight)
{
    auto cur_edge = std::make_shared<Edge>(vertex_list[node_a], vertex_list[node_b], weight);
    vertex_list[node_a]->add_edge(cur_edge);
    vertex_list[node_b]->add_edge(cur_edge);
    edge_list.push_back(std::move(cur_edge));
}


const Path Graph::find_path(int start, int end) const
{
    int cur_i = start;
    double cur_weight = 0.0;
    std::unordered_set<int> closed_set = {};  // only need to remember which nodes have been visited
    MinHeap<int> open_weight;  // Min heap to make finding minimum value fast
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

    void run_week3_homework()
    {
        double orig_density = density;
        density = 0.4;
        test_shortest_path();
        density = 0.2;
        test_shortest_path();
        density = orig_density;
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
    enum class TestTypes {TestPath, Week3Homework, TestRead};

    TestTypes test_to_run = TestTypes::TestRead;
    if (test_to_run == TestTypes::Week3Homework) {
        setup.run_week3_homework();
    }
    else if (test_to_run == TestTypes::Week3Homework) {
        // Run debugging setup
        setup.test_shortest_path();
    }
    else if (test_to_run == TestTypes::TestRead) {
        std::string test_fname("./sample_graph.txt");
        std::cout << "Attempting to read \"" << test_fname << "\"" << std::endl;
        Graph read_test_graph(test_fname);
        std::cout << read_test_graph;
    }
    else {
        std::cout << "Unknown Run Setup" << std::endl;
    }
    return 0;
}
