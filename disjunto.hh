#include "nodeedge.hh"
#include <iostream>
#include <unordered_map>
using namespace std;

// template <typename T, typename N>
class disjoint {
private:
  unordered_map<string, string> Parent;
  unordered_map<string, int> Rank;

public:
  /*disjoint() {

    char universe[] = {'a', 'b', 'c', 'd', 'e'};
    for (char x : universe) {
      Parent[x] = x;
    }
    Parent['d'] = 'b';
  }*/

  string find(string vertex) {
    if (Parent[vertex] == vertex)
      return vertex;
    else
      return find(Parent[vertex]);
  }

  // void Union(Node *root1, Node *root2) { Parent[root1] = root2; }
  void UnionbyRank(string x, string y) {
    string xRoot = find(x); // Obtengo la raiz de la componente del vértice X
    string yRoot = find(y); // Obtengo la raiz de la componente del vértice Y
    if (Rank[xRoot] > Rank[yRoot]) { // en este caso la altura de la componente
      // del vértice X es mayor que la altura de
      // la componente del vértice Y.
      Parent[yRoot] = xRoot;
      // el padre de ambas componentes será el de mayor altura
    } else { // en este caso la altura de la componente del vértice Y es mayor o
      // igual que la de X
      Parent[xRoot] = yRoot;
      // el padre de ambas componentes será el de mayor altura
      if (Rank[xRoot] == Rank[yRoot]) { // si poseen la misma altura
        Rank[yRoot]++;                  // incremento el rango de la nueva raíz
      }
    }
  }

  void MakeSet(Node *vertex) {
    Parent[vertex->getId()] = vertex->getId();
    Rank[vertex->getId()] = 0;
  }

  void print() {
    cout << "My hash table contains: " << endl;
    for (auto it = Parent.begin(); it != Parent.end(); ++it) {
      cout << it->first << ": " << it->second << endl;
    }
  }
};
