/**
 * Graph Class for Dijkstra's Algorithm
 */


#ifndef _GRAPH_H
#define _GRAPH_H


#include <iostream>
#include <memory>
#include <vector>
#include <list>
#include <queue>
#include <unordered_map>
#include <unordered_set>


namespace dijkstra {

    /**
     * Queue compare function to create min-heap
    */
    template <class T, typename HASH=std::hash<T>>
    bool compare_second(std::pair<T, double> &a, std::pair<T, double> &b) {
        return a.second > b.second;
    }

    /**
     * Min Heap Implementation
     * 
     * This uses a priority queue to keep track of the minimum weight and a dict to make accessing the weight easy.
    */
    template <class T, typename HASH=std::hash<T>>
    class MinHeap {
        public:

            /**
             * Helper queue type to make MinHeap more readable.  Data is <index, weight>.
            */
            typedef std::pair<T, double> queue_val;

            /**
             * Constructor for empty min heap.
            */
            MinHeap();

            /**
             * Set index and value.
            */
            void set(T index, double value);

            /**
             * Get value associated with index.
            */
            double get(T index) const;

            /**
             * Accessor-only indexing (can't set)
            */
            double operator[](T index) const;

            /**
             * Helper function to see if the index is in the min-heap
            */
            double contains(T index) const;

            /**
             * Remove index from heap
            */
            void erase(T index);

            /**
             * Pop the top of the min-heap
            */
            queue_val get_min();

            /**
             * Get sizre of min-heap
            */
            int size() const;

            friend std::ostream &operator<<(std::ostream &os, MinHeap &min_heap)
            {
                os << "MinHeap(";
                for (auto item: min_heap.val_map) {
                    os << "<" << item.first << ", " << item.second << ">, ";
                }
                os << ")";
                return os;
            }

        private:

            std::priority_queue<queue_val, std::vector<queue_val>, decltype(&compare_second<T>)> queue;
            std::unordered_map<T, double, HASH> val_map;

    };


    // Circular referencing shenanigans.  Let the compiler know Edge will be defined now.
    class Edge;


    /**
     * Vertex class that holds Edge objects connecting to other Vertices
    */
    class Vertex {
        public:

            /**
             * Create vertex indexed by `id`
            */
            Vertex(int id);

            /**
             * Get the ID of this vertex
            */
            int get_id() const;

            /**
             * Add an edge to this vertex
            */
            void add_edge(std::shared_ptr<Edge> edge);

            /**
             * Check to see if this vertex is connected to another vertex
            */
           bool is_connected_to(int id);

            /**
             * Get a list of vertices that are connected and their path weights
            */
            const std::vector<std::pair<int, double>> get_edges() const;

            /**
             * Get the "A" side of all edges attached to this node.  Useful if A vs B is meaningful (e.g. in vs out)
            */
            const std::vector<std::pair<int, double>> get_edges_a() const;

            /**
             * Get the "B" side of all edges attached to this node.  Useful if A vs B is meaningful (e.g. in vs out)
            */
            const std::vector<std::pair<int, double>> get_edges_b() const;

            friend std::ostream &operator<<(std::ostream &os, Vertex const &vert);
        private:
            int id;
            std::vector<std::shared_ptr<Edge>> conn_list;
    };


    /**
     * Edge class that connects two vertices
    */
    class Edge {
        public:

            /**
             * Create edge from Vertex A to Vertex B, with associated weight
            */
            Edge(std::shared_ptr<Vertex> vertex_a, std::shared_ptr<Vertex> vertex_b, double weight);

            /**
             * Get pointer to vertex A
            */
            const Vertex &get_a() const;

            /**
             * Get pointer to vertex B
            */
            const Vertex &get_b() const;

            /**
             * "Travel" this edge from input vertex, returning the vertex on the other side of the edge.
            */
            const Vertex &get_other(const Vertex &thing) const;

            /**
             * Get edge weight
            */
            double get_weight() const;

            friend std::ostream &operator<<(std::ostream &os, Edge const &edge)
            {
                os << "Edge: " << edge.vertex_a->get_id() << " <-> " << edge.vertex_b->get_id()
                    << ", weight" << edge.weight;
                return os;
            }
        private:
            std::shared_ptr<Vertex> vertex_a;
            std::shared_ptr<Vertex> vertex_b;
            double weight;
    };


    // Need to define this below Edge, since the compiler doesn't know what Edge is when Vertix is declared
    std::ostream &operator<<(std::ostream &os, Vertex const &vert)
    {
        os << "Vertex " << vert.id << "\n";
        for (auto edge: vert.conn_list) {
            os << " - " << edge->get_other(vert).get_id() << " (weight=" << edge->get_weight() << ")\n";
        }
        return os;
    }


    /**
     * Output path class that holds the path weight and node index list.
    */
    class Path {
        public:

            /**
             * Create output path with total weight
            */
            Path(double weight);

            /**
             * Get weight of path
            */
            double get_weight() const;

            /**
             * Push node into front of path
            */
            void push_front(int id);

            /**
             * Push node into back of path
            */
            void push_back(int id);

            /**
             * Get the path
            */
            const std::list<int> get_path() const;

            friend std::ostream &operator<<(std::ostream &os, Path const &path)
            {
                os << "Path: weight=" << path.weight << "\n - ";
                for (auto elem_i : path.path) {
                    os << elem_i << " -> ";
                }
                os << std::endl;
                return os;
            }
        private:
            double weight;
            std::list<int> path;
    };


    class Tree {
        public:

            /**
             * Create empty tree
             */
            Tree();

            /**
             * Create tree with specified head
            */
            Tree(int head_id);

            /**
             * Set head node
             */
            void set_head(int head_id);

            /**
             * Get head node
            */
            const Vertex &get_head() const;

            /**
             * Add child to tree
            */
            void add_child(int parent_id, int child_id, double weight);

            /**
             * Get total weight
            */
            double get_total_weight() const;

            friend std::ostream &operator<<(std::ostream &os, Tree const &tree)
            {
                os << "Tree: w=weight, n=node\n";
                tree.dump_children(os, tree.head->get_id(), 0.0, 0);
                os << "Total Weight = " << tree.get_total_weight() << std::endl;
                return os;
            }

        private:
            std::shared_ptr<Vertex> head;
            std::unordered_map<int, std::shared_ptr<Vertex>> vertex_map;
            std::vector<std::shared_ptr<Edge>> edge_list;

            /**
             * Helper function for printing tree with tabs
            */
            void dump_children(std::ostream &os, int node_i, double weight, int n_indent) const;

    };


    /**
     * Graph path for dijkstra's algorithm
    */
    class Graph {
        public:

            /**
             * Create empty graph with n_nodes and no edges
            */
            Graph(int n_nodes);

            /**
             * Create random graph for testing
             * 
             * @param[in] n_nodes  Number of nodes in path
             * @param[in] density  Proportion of connected nodes
             * @param[in] dist_min  Minimum distance between nodes (uniform)
             * @param[in] dist_max  Maximum distance between nodes (uniform)
            */
            Graph(int n_nodes, double density, double dist_min, double dist_max);

            /**
             * Create graph from file.
             * 
             * First line of file shall be the number of nodes.
             * Each other line is an edge and has three numbers (node_i, node_j, cost)
            */
            Graph(std::string graph_fname);

            /**
             * Helper function to add edge to graph
            */
            void add_edge(int node_a, int node_b, double weight);

            /**
             * Find path between node indices using dijkstra's algorithm
            */
            const Path find_path(int start, int end) const;

            /**
             * Generate minimum spanning tree, starting at specified node
            */
            const Tree find_mst_prim(int head) const;
            
            friend std::ostream &operator<<(std::ostream &os, Graph const &graph)
            {
                os << "Graph: n_nodes=" << graph.n_nodes << "\n";
                os << "^^^^^^^^^^^^^^^^^^^^\n";
                for (auto vert : graph.vertex_list) {
                    os << *vert;
                }
                os << "vvvvvvvvvvvvvvvvvvvv" << std::endl;
                return os;
            }

        private:
            int n_nodes;
            std::vector<std::shared_ptr<Vertex>> vertex_list;
            std::vector<std::shared_ptr<Edge>> edge_list;

            /**
             * Helper function to create node list
            */
            void create_node_list(int n_nodes);

    };

}

#endif
