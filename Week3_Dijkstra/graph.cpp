/* Graph Class for Dijkstra's Algorithm
 * 
 * Week 3 assignment - C++ for C programmers
 * 
 */


#include <random>
#include <algorithm>
#include "graph.h"


using namespace dijkstra;


Vertex::Vertex(int id): id(id), conn_list(0)
{}


int Vertex::get_id() const
{
    return id;
}


void Vertex::add_edge(std::shared_ptr<Edge> edge)
{
    conn_list.push_back(edge);
}


std::unordered_map<int, double> Vertex::get_edges() const
{
    std::unordered_map<int, double> out_map;
    for (auto edge : conn_list) {
        int other_i = edge->get_other(*this)->get_id();
        double weight = edge->get_weight();
        out_map[other_i] = weight;
    }
    return out_map;
}


Edge::Edge(Vertex &vertex_a, Vertex &vertex_b, double weight):
    vertex_a(&vertex_a), vertex_b(&vertex_b), weight(weight)
{}


std::shared_ptr<Vertex> Edge::get_a() const
{
    return vertex_a;
}


std::shared_ptr<Vertex> Edge::get_b() const
{
    return vertex_b;
}


std::shared_ptr<Vertex> Edge::get_other(const Vertex &thing) const
{
    if(thing.get_id() == vertex_a->get_id()) {
        return vertex_a;
    }
    else if (thing.get_id() == vertex_b->get_id()) {
        return vertex_b;
    }
    else {
        return nullptr;
    }
}


double Edge::get_weight() const
{
    return weight;
}


Path::Path(double weight=0): weight(0), path()
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
    std::unordered_map<int, double> node_weight = {};
    std::unordered_map<int, int> prev_dict = {};

    node_weight[start] = 0.0;

    // Limit the number of nodes we check to the number of nodes, in case I do something
    // silly
    for (int iter_i=0; iter_i < vertex_list.size(); ++iter_i) {
        // We're done if we're at the end
        if (cur_i == end) {
            break;
        }
        // Get the current weight and the edges out of this vertex
        double cur_weight = node_weight[cur_i];
        const std::unordered_map<int, double> out_edges = vertex_list[cur_i]->get_edges();
        // Check each edge leaving the vertex
        for (auto edge: out_edges) {
            // Don't bother checking if we've already closed the other side of the edge
            if(closed_set.find(edge.first) != closed_set.end()) {
                continue;
            }
            double new_weight = cur_weight + edge.second;
            if (node_weight.find(edge.first) != node_weight.end()) {
                // vertex has already been visited, update if this has a lower cost than previous visit.
                if(new_weight < node_weight[edge.first]) {
                    node_weight[edge.first] = new_weight;
                    prev_dict[edge.first] = cur_i;
                }
            }
            else {
                // Never visited this vertex, so update
                node_weight[edge.first] = new_weight;
                prev_dict[edge.first] = cur_i;
            }
        }
        // Done with current vertex, add to closed set and remove from possible nodes.
        closed_set.insert(cur_i);
        node_weight.erase(cur_i);
        // Find next vertex to work from (lowest weight thus far)
        std::pair<int, double> cur_min_pair = {-1, INFINITY};
        for (auto test_weight_pair: node_weight) {
            // Checking closed should be unnecessary?
            // if (closed_set.find(test_weight_pair.first) != closed_set.end()) {
            //     continue;
            // }
            if (test_weight_pair.second < cur_min_pair.second) {
                cur_min_pair = test_weight_pair;
            }
        }
        // No more nodes to check, so end must not be reachable from start.
        if (cur_i == -1) {
            // Not sure how to do exceptions yet, so I'm just going to return an empty path
            return Path(-1);
        }
        cur_i = cur_min_pair.first;
    }
    // If we're here, we found a path
    Path out_path(node_weight[end]);
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
