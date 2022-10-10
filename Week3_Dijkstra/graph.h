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

    typedef std::pair<int, double> queue_val;

    bool compare_second(queue_val &a, queue_val &b) {
        return a.second > b.second;
    }

    class MinHeap {
        public:
            MinHeap();
            void set(int index, double value);
            double get(int index) const;
            double operator[](int index) const;
            double contains(int index) const;
            void erase(int index);
            queue_val get_min();
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


    class Vertex {
        public:
            Vertex(int id);
            int get_id() const;
            void add_edge(std::shared_ptr<Edge> edge);
            const std::vector<std::pair<int, double>> get_edges() const;
            friend std::ostream &operator<<(std::ostream &os, Vertex const &vert);
        private:
            int id;
            std::vector<std::shared_ptr<Edge>> conn_list;
    };


    class Edge {
        public:
            Edge(std::shared_ptr<Vertex> vertex_a, std::shared_ptr<Vertex> vertex_b, double weight);
            const Vertex &get_a() const;
            const Vertex &get_b() const;
            const Vertex &get_other(const Vertex &thing) const;
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


    class Path {
        public:
            Path(double weight);
            double get_weight() const;
            void push_front(int id);
            void push_back(int id);
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


    class Graph {
        public:
            Graph(int n_nodes, double density, double dist_min, double dist_max);
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
