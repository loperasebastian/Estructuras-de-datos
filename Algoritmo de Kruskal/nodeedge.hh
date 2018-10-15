//@brief Definition of a particular node in the graph.
#include <iostream>

using namespace std;
class Node {
private:
  string id_;

public:
  Node(void) {}
  Node(string id) { id_ = id; }
  string getId(void) const { return id_; }
  void setId(string id) { id_ = id; }
};

/**
 * @brief Definition of a particular edge in the graph.
 */
class Edge {
private:
  Node *start_;
  Node *end_;
  int weight_;

public:
  Edge(void) { weight_ = 0; }

  void setStart(Node *s) { start_ = s; }
  void setEnd(Node *e) { end_ = e; }

  Node *getStart(void) const { return start_; }
  Node *getEnd(void) const { return end_; }

  void setWeight(int w) { weight_ = w; }
  int getWeight(void) const { return weight_; }
};
