/* Heap Test
 */


#include <iostream>
#include <random>
#include <vector>
#include <queue>
#include <memory>
#include <unordered_map>
#include <numeric>


typedef std::pair<int, double> queue_val;


bool compare_second(queue_val &a, queue_val &b) {
    return a.second > b.second;
}


class MinHeap {
    public:
        MinHeap(): queue(compare_second), val_map() {}
        void set(int index, double value)
        {
            val_map[index] = value;
            queue.push(queue_val(index, value));
        }
        double get(int index) const
        {
            // Note unordered_map operator[] is not const because it will add an element if it is not found.
            return val_map.at(index);
        }
        double operator[](int index) const
        {
            return get(index);
        }
        queue_val get_min() {
            queue_val out;
            while (!queue.empty()) {
                out = queue.top();
                queue.pop();
                if (val_map.find(out.first) != val_map.end() && val_map[out.first] == out.second) {
                    val_map.erase(out.first);
                    return out;
                }
            }
            out = {-1, 0.0};
            return out;
        }
        int size() const
        {
            return val_map.size();
        }
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


int main()
{
    // Initialize random number generator for the graph
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> edge_dist(0.0, 1.0);

    MinHeap heap;
    std::cout << "Creating Heap\n";
    for (int val_i=0; val_i < 6; ++val_i) {
        double val = edge_dist(generator);
        std::cout << " - <" << val_i << ", " << val << ">\n";
        heap.set(val_i, val);
    }
    std::cout << std::endl;

    std::cout << "Popping Output\n";
    for (int val_i=0; val_i < 3; ++val_i) {
        queue_val out = heap.get_min();
        std::cout << " - <" << out.first << ", " << out.second << ">\n";
    }
    std::cout << std::endl;

    for (int val_i=0; val_i < 6; ++val_i) {
        double val = edge_dist(generator);
        std::cout << " - <" << val_i << ", " << val << ">\n";
        heap.set(val_i, val);
    }
    std::cout << std::endl;

    std::cout << "Popping Output\n";
    for (int val_i=0; val_i < 6; ++val_i) {
        queue_val out = heap.get_min();
        std::cout << " - <" << out.first << ", " << out.second << ">\n";
    }
    std::cout << std::endl;

    return 0;
}
