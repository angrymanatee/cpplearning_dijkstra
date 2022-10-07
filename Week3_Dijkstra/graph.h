/**
 * Graph Class for Dijkstra's Algorithm
 */


#ifndef _GRAPH_H
#define _GRAPH_H


#include <memory>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>


namespace dijkstra {

    class Vertex {
        public:
            Vertex(int id);
            int get_id() const;
            void add_edge(std::shared_ptr<Edge> edge);
            const std::vector<std::pair<int, double>> get_edges() const;
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
        private:
            std::shared_ptr<Vertex> vertex_a;
            std::shared_ptr<Vertex> vertex_b;
            double weight;
    };

    class Path {
        public:
            Path(double weight=0);
            double get_weight() const;
            void push_front(int id);
            void push_back(int id);
            const std::list<int> get_path() const;
        private:
            double weight;
            std::list<int> path;
    };

    class Graph {
        public:
            Graph(int n_nodes, double density, double dist_min, double dist_max);
            const Path find_path(int start, int end) const;
        private:
            int n_nodes;
            double density;
            double dist_min;
            double dist_max;
            std::vector<std::shared_ptr<Vertex>> vertex_list;
    };

}

#endif
