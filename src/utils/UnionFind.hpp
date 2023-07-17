#pragma once
#include <algorithm>
#include <assert.h>
#include <vector>
struct UnionFind {
  std::vector<int> parent;
  std::vector<int> size;
  UnionFind(int s) : parent(s), size(s) {
    for (int &p : parent) {
      parent[p] = p;
      size[p] = 1;
    }
  }
  int find_set(int v) {
    if (v == parent[v])
      return v;
    return parent[v] = find_set(parent[v]);
  }

  void union_sets(int a, int b) {
    a = find_set(a);
    b = find_set(b);
    if (a != b) {
      if (size[a] < size[b])
        std::swap(a, b);
      parent[b] = a;
      size[a] += size[b];
    }
  }
};
