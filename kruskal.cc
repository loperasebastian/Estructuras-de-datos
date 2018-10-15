#include "disjunto.hh"
#include "graph.hh"
#include <iostream>
#include <vector>

using namespace std;

void Kruskal(Graph<Node, Edge> *g, vector<Edge *> &a, vector<Node *> &b) {
  vector<Edge *> MST;
  disjoint ds;
  int peso = 0;
  for (Node *d : b) {
    ds.MakeSet(d);
  }

  // ds.print();
  for (Edge *e : a) {
    if (ds.find(e->getStart()->getId()) != ds.find(e->getEnd()->getId())) {
      MST.push_back(e);
      ds.UnionbyRank(e->getStart()->getId(), e->getEnd()->getId());
    }
  }

  for (int i = 0; i < MST.size(); i++) {

    cout << "arco " << i + 1 << "   " << MST[i]->getStart()->getId()
         << MST[i]->getEnd()->getId() << MST[i]->getWeight() << endl;
    peso = peso + MST[i]->getWeight();
  }
  cout << peso << endl;
}

int main() {
  Graph<Node, Edge> *g = new Graph<Node, Edge>();

  vector<Edge *> a;
  vector<Node *> b;
  g->addNode("A");
  g->addNode("B");
  g->addNode("C");
  g->addNode("D");
  g->addNode("E");
  g->addNode("F");
  g->addNode("G");
  g->addNode("H");
  g->addNode("I");
  g->addUndirectedEdge("A", "B", 4);
  g->addUndirectedEdge("A", "H", 8);
  g->addUndirectedEdge("B", "C", 8);
  g->addUndirectedEdge("B", "H", 11);
  g->addUndirectedEdge("H", "G", 1);
  g->addUndirectedEdge("H", "I", 7);
  g->addUndirectedEdge("I", "C", 2);
  g->addUndirectedEdge("I", "G", 6);
  g->addUndirectedEdge("C", "F", 4);
  g->addUndirectedEdge("C", "D", 7);
  g->addUndirectedEdge("G", "F", 2);
  g->addUndirectedEdge("F", "E", 10);
  g->addUndirectedEdge("D", "E", 9);
  g->addUndirectedEdge("D", "F", 14);

  for (Edge *b : g->getEdges()) {
    a.push_back(b);
  }
  for (Node *c : g->getNodes()) {
    b.push_back(c);
  }
  /*for (int i = 0; i < b.size(); i++) {
    cout << "Nodo" << i + 1 << "  " << b[i]->getId() << endl;
  }
  for (int i = 0; i < a.size(); i++) {

    cout << "arco " << i + 1 << "   " << a[i]->getStart()->getId()
         << a[i]->getEnd()->getId() << a[i]->getWeight() << endl;
  }
*/
  g->quicksort(a, 0, a.size() - 1);
  /*for (int i = 0; i < a.size(); i++) {

    cout << a[i]->getWeight() << a[i]->getStart()->getId()
         << a[i]->getEnd()->getId() << endl;
  }*/

  Kruskal(g, a, b);
  return 0;
}
