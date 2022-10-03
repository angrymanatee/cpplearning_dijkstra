/* Graph Class for Dijkstra's Algorithm
 * 
 * Week 3 assignment - C++ for C programmers
 * 
 */


#include <vector>


using namespace dijkstra;

class DijEdge {
    public:
        DijEdge(DijVertex &from, DijVertex &to, double weight):
            from(from), to(to), weight(weight) {};
        DijVertex &from;
        DijVertex &to;
        double weight;
};

class DijVertex {
    public:
        DijVertex(int id): id(id), conn_list() {};
        int id;
        std::vector<DijEdge> conn_list;
};
