//  Convert this program to C++
//  change to C++ io
//  change to one line comments
//  change defines of constants to const
//  change array to vector<>
//  inline any short function


#include <vector>
#include <iostream>


const int N = 40;
using namespace std;


// Sum the first n elements in a vector
template <class T> inline void sum(T &p, int n, const vector<T> &d) {
  p = 0;
  for(int i = 0; i < n; ++i) {
    p = p + d[i];
  }
}


int main() {
  int accum = 0;
  vector<int> data(N, 0);  // Initialize to 0 since it won't effect summing
  for(int i = 0; i < N; ++i) {
    data[i] = i;
  }
  sum(accum, N, data);
  cout << "sum is " << accum << endl;
  return 0;
}
