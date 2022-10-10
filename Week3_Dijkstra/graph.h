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
     * Helper queue type to make MinHeap more readable.  Data is <index, weight>.
    */
    typedef std::pair<int, double> queue_val;

    /**
     * Queue compare function to create min-heap
    */
    bool compare_second(queue_val &a, queue_val &b) {
        return a.second > b.second;
    }

    /**
     * Min Heap Implementation
     * 
     * This uses a priority queue to keep track of the minimum weight and a dict to make accessing the weight easy.
    */
    class MinHeap {
        public:

            /**
             * Constructor for empty min heap.
            */
            MinHeap();

            /**
             * Set index and value.
            */
            void set(int index, double value);

            /**
             * Get value associated with index.
            */
            double get(int index) const;

            /**
             * Accessor-only indexing (can't set)
            */
            double operator[](int index) const;

            /**
             * Helper function to see if the index is in the min-heap
            */
            double contains(int index) const;

            /**
             * Remove index from heap
            */
            void erase(int index);

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
            std::priority_queue<queue_val, std::vector<queue_val>, decltype(&compare_second)> queue;
            std::unordered_map<int, double> val_map;
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
             * Get a list of vertices that are connected and their path weights
            */
            const std::vector<std::pair<int, double>> get_edges() const;

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


    /**
     * Graph path for dijkstra's algorithm
    */
    class Graph {
        public:
        
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
             * Find path between node indices using dijkstra's algorithm
            */
            const Path find_path(int start, int end) const;
            
            friend std::ostream &operator<<(std::ostream &os, Graph const &graph)
            {
                os << "Graph: n_nodes=" << graph.n_nodes << ", density=" << graph.density
                    << ", dist=(" << graph.dist_min << ", " << graph.dist_max << ")\n";
                for (auto vert : graph.vertex_list) {
                    os << *vert;
                }
                os << std::endl;
                return os;
            }
        private:
            int n_nodes;
            double density;
            double dist_min;
            double dist_max;
            std::vector<std::shared_ptr<Vertex>> vertex_list;
    };

}

#endif
