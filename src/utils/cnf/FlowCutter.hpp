#pragma once

#include "cell.h"
#include "id_multi_func.h"
#include "src/specs/cnf/SpecManagerCnf.hpp"

struct TreeDecomposition {
  std::vector<std::vector<d4::Var>> adj;
  std::vector<std::vector<int>> bags;
  std::vector<int> cnt;
  int width;

  int find_centroid() {
    cnt.clear();
    cnt.resize(bags.size());
    for (int i = 0; i < bags.size(); i++) {
      for (int j = 0; j < adj[i].size(); j++) {
        cnt[j] += 1;
      }
    }
  }
};

struct FlowCutter {
  ArrayIDIDFunc head, tail;
  ArrayIDIDFunc preorder, inv_preorder;
  int best_bag_size = std::numeric_limits<int>::max();
  long long timeout = 30000;
  int seed = 42;
  TreeDecomposition best_tree;
  d4::SpecManagerCnf *om;
  FlowCutter(d4::SpecManagerCnf *specs);
  void update_best_order(ArrayIDIDFunc &tail, ArrayIDIDFunc &head,
                         const ArrayIDIDFunc &order);
  void update_best_partition(ArrayIDIDFunc &tail, ArrayIDIDFunc &head,
                             const ArrayIDIDFunc &to_input_node_id,
                             const std::vector<Cell> &cell_list);
  void decompose(std::vector<d4::Var> &vars);
  void test_new_order(const ArrayIDIDFunc &order);
  void create_primal(std::vector<d4::Var> vars);
  int compute_max_bag_size_of_order(const ArrayIDIDFunc &order);
};
