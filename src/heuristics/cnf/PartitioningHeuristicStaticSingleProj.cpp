/*
 * d4
 * Copyright (C) 2020  Univ. Artois & CNRS
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "PartitioningHeuristicStaticSingleProj.hpp"

#include <ostream>

#include "src/utils/AtMost1Extractor.hpp"
#include "src/utils/EquivExtractor.hpp"
#include "src/utils/Proj.hpp"
#include <numeric>
#define LOG 0

namespace d4 {

namespace hyper_util {

void saveHyperGraph(std::vector<std::vector<unsigned>> &savedHyperGraph,
                    std::vector<int> &savedCost, HyperGraph &h) {
  for (auto edge : h) {
    savedHyperGraph.push_back(std::vector<unsigned>());
    std::vector<unsigned> &tmp = savedHyperGraph.back();
    for (auto v : edge)
      tmp.push_back(v);
    savedCost.push_back(h.getCost()[edge.getId()]);
  }
}
void setHyperGraph(std::vector<std::vector<unsigned>> &savedHyperGraph,
                   std::vector<int> &savedCost, std::vector<unsigned> &indices,
                   HyperGraph &hypergraph) {

  unsigned *edges = hypergraph.getEdges();
  hypergraph.setSize(0);

  for (auto idxEdge : indices) {
    std::vector<unsigned> &tmp = savedHyperGraph[idxEdge];
    if (!tmp.size())
      continue;

    *edges = tmp.size();
    for (unsigned i = 0; i < tmp.size(); i++)
      edges[i + 1] = tmp[i];
    edges += *edges + 1;
    hypergraph.getCost()[hypergraph.getSize()] = savedCost[idxEdge];

    // hypergraph.getCost()[hypergraph.getSize()] =
    //     savedCost[idxEdge] > 1 ? std::min<int>(2, indices.size()*0.1) : 1;
    hypergraph.incSize();
  }
}

} // namespace hyper_util

/**
   Constructor.

   @param[in] vm, the option list.
   @param[in] s, a wrapper on a solver.
   @param[in] om, a structure manager.
*/
PartitioningHeuristicStaticSingleProj::PartitioningHeuristicStaticSingleProj(
    po::variables_map &vm, WrapperSolver &s, SpecManager &om, std::ostream &out)
    : PartitioningHeuristicStaticSingleProj(
          vm, s, om, dynamic_cast<SpecManagerCnf &>(om).getNbClause(),
          dynamic_cast<SpecManagerCnf &>(om).getNbVariable(),
          dynamic_cast<SpecManagerCnf &>(om).getSumSizeClauses(), out) {

} // constructor

/**
   Constructor.

   @param[in] vm, the option list.
   @param[in] s, a wrapper on a solver.
   @param[in] om, a structure manager.
   @param[in] nbClause, the number of clauses.
   @param[in] nbVar, the number of variables.
   @param[in] sumSize, which give the number of literals.
 */
PartitioningHeuristicStaticSingleProj::PartitioningHeuristicStaticSingleProj(
    po::variables_map &vm, WrapperSolver &s, SpecManager &om, int nbClause,
    int nbVar, int sumSize, std::ostream &out)
    : PartitioningHeuristicStaticSingle(vm, s, om, nbClause, nbVar, sumSize,
                                        out) {

  m_score = ScoringMethod::makeScoringMethod(vm, om, s, out);
  m_max_cut_ratio = vm["partitioning-heuristic-max-cut-ratio"].as<float>();
} // constructor

/*
std::vector<int> get_cut(HyperGraph &h, const std::vector<int> &part,
                         const std::vector<int> &map) {
  std::vector<int> cut;
  for (auto &e : h) {
    int p = map[part[*e.begin()]];
    for (auto v : e) {
      if (p != map[part[v]]) {
        cut.push_back(e.getId());
        break;
      }
    }
  }
  return cut;
}
void reduce_multi_cut(HyperGraph &h,std::vector<bool> removed, std::vector<int>
&part, int part_cnt) { size_t min_cost = std::numeric_limits<size_t>::max();
  std::vector<int> map(part_cnt / 2, 0), min_map;
  map.resize(part_cnt, 1);
  while (std::next_permutation(map.begin(), map.end())) {
    auto cut = get_cut(h, part, map);
    size_t cost = 0;
    for (auto e : cut) {
      cost += h.getCost()[e];
    }
    if (cost < min_cost) {
      min_cost = cost;
      min_map = map;
    }
  }
  for (auto i : get_cut(h, part, min_map)) {
  }
  std::cout<<"\n";
  for (auto &p : part) {
    p = map[p];
  }
}
*/
/**
   Split and assign variables.

   @param[in] indicesFirst, the first parition.
   @param[in] indicesSecond, the second partition.
   @param[in] mappingVar, to get the variable associate with the index.
   @param[in] cutIsempty, specify if the partition that generates the two set of
   indices come from an empty cut set.
   @param[out] stack, the current stack of set of variables (will receive
   indicesFirst and indicesSecond if their size is large enough).
   @param[out] level, the current level where are assigned the variables in
   their bucket.
*/

/**
   Search a decomposition tree regarding a component.

   @param[in] component, the set of varaibles the problem is constructed on.
   @param[in] equivClass, the equivalence class for each variable.
   @param[in] equivVar, the list of equivalences.
   @param[out] bucketNumber, the decomposition tree in term of index.
*/

std::vector<unsigned> get_cut(HyperGraph &h, const std::vector<int> &part,
                              const std::vector<int> &map,
                              const std::vector<int> &base_map) {

  std::vector<unsigned> cut;
  for (auto &e : h) {
    int base = base_map[part[*e.begin()]];
    if (base == -1) {
      continue;
    }
    int p = map[base];
    bool in_cut = false;
    for (auto v : e) {
      int base = base_map[part[v]];
      if (base == -1) {
        in_cut = false;
        break;
      }
      if (p != map[base]) {
        in_cut = true;
      }
    }
    if (in_cut) {
      cut.push_back(e.getId());
    }
  }
  return cut;
}

void split(HyperGraph &h, const std::vector<int> &part,
           const std::vector<int> &map, const std::vector<int> &base_map,
           std::vector<unsigned> &cut, std::vector<unsigned> &indices0,
           std::vector<unsigned> &indices1) {
  for (auto &e : h) {
    int base = base_map[part[*e.begin()]];
    if (base == -1) {
      continue;
    }
    int p = map[base];
    bool in_cut = false;
    for (auto v : e) {
      int base = base_map[part[v]];
      if (base == -1) {
        in_cut = false;
        break;
      }
      if (p != map[base]) {
        in_cut = true;
      }
    }
    if (in_cut) {
      cut.push_back(e.getId());
    } else if (p == 0) {
      indices0.push_back(e.getId());
    } else {
      indices1.push_back(e.getId());
    }
  }
}

std::vector<int> optimize_multi_cut(HyperGraph &h, const std::vector<int> &part,
                                    const std::vector<int> &base_map,
                                    const int part_cnt) {
  size_t min_cost = std::numeric_limits<size_t>::max();
  std::vector<int> map(part_cnt / 2, 0), min_map;
  map.resize(part_cnt, 1);
  while (std::next_permutation(map.begin(), map.end())) {
    auto cut = get_cut(h, part, map, base_map);
    size_t cost = 0;
    for (auto e : cut) {
      cost += h.getCost()[e];
    }
    if (cost < min_cost) {
      min_cost = cost;
      min_map = map;
    }
  }
  return map;
}
std::vector<int> reduce_map(std::vector<int> base_map,
                            std::vector<int> reduction, int x) {
  std::vector<int> out(base_map.size());
  int z = 0;
  for (int k = 0; k < out.size(); k++) {
    if (base_map[k] == -1) {
      out[k] = -1;
    } else {
      int p = reduction[base_map[k]];
      if (p == x) {
        out[k] = z;
        z++;
      } else {
        out[k] = -1;
      }
    }
  }
  return out;
}

void PartitioningHeuristicStaticSingleProj::handle_multicut(
    HyperGraph &hypergraph, std::vector<std::vector<unsigned>> &hypergraph_list,
    std::vector<int> &partition, std::vector<unsigned> &mappingEdge,
    std::vector<Var> &mappingVar, std::vector<Strata> &stack, unsigned &level) {
  struct MultiCutLvl {
    std::vector<int> part_map;
    int depth;
    unsigned fatherId;
  };
  std::vector<MultiCutLvl> lvls;
  std::vector<int> inital_map(m_multi_lvl_cuts);
  std::iota(inital_map.begin(), inital_map.end(), 0);
  lvls.push_back({inital_map, 0, stack.back().fatherId});
  auto push_new_lvl = [&](std::vector<unsigned> &indices,
                          std::vector<int> base_map, std::vector<int> reduction,
                          int depth, unsigned fatherId, int x) {
    if (indices.size() > LIMIT) {
      if ((1 << depth) == m_multi_lvl_cuts) {
        stack.push_back(Strata{fatherId, indices});

      } else {
        lvls.push_back(
            MultiCutLvl{.part_map = reduce_map(base_map, reduction, x),
                        .depth = depth,
                        .fatherId = fatherId});
      }
    } else {
      assignLevel(hypergraph_list, fatherId, indices, mappingVar, level);
    }
  };
  while (lvls.size()) {
    auto lvl = std::move(lvls.back());
    auto fatherId = lvl.fatherId;
    lvls.pop_back();
    std::cout << "Map: ";
    for (auto m : lvl.part_map) {
      std::cout << m << " ";
    }
    std::cout << std::endl;
    std::vector<int> reduction_map =
        optimize_multi_cut(hypergraph, partition, lvl.part_map,
                           m_multi_lvl_cuts / (1 << (lvl.depth)));
    std::vector<unsigned> cut, indices0, indices1;
    split(hypergraph, partition, reduction_map, lvl.part_map, cut, indices0,
          indices1);
    unsigned currentId = (cut.size()) ? level : fatherId;
    if (cut.size()) {
#if LOG
      m_log << "Lvl: " << level << ": ";
      int c = 0;
      for (auto i : cut) {
        Var v = m_hypergraph.getCost()[i];
        c += v;
        m_log << v << ":";
        m_log << mappingVar[i] << " ";
      }
      m_log << "cost: " << c << "\n";
#endif
      setBucketLevelFromEdges(hypergraph_list, cut, mappingVar, level);
      assert(fatherId < m_levelInfo.size());
      m_levelInfo[fatherId].separatorLevel = level;
      level++;
      m_levelInfo.push_back({level, (unsigned)cut.size()});
    } else {
      // special case 1.
      if (!indices0.size() && indices1.size()) {
        assignLevel(hypergraph_list, currentId, indices1, mappingVar, level);
        continue;
      }
      // special case 2.
      if (!indices1.size() && indices0.size()) {
        assignLevel(hypergraph_list, currentId, indices0, mappingVar, level);
        continue;
      }
    }
    push_new_lvl(indices1, lvl.part_map, reduction_map, lvl.depth + 1,
                 currentId, 1);
    push_new_lvl(indices0, lvl.part_map, reduction_map, lvl.depth + 1,
                 currentId, 0);
  }
}

void probe(WrapperSolver &s, SpecManager &specs, ScoringMethod *score,
           ProjVars &vars, std::vector<float> &units_cover,
           std::vector<bool> &cut, int depth) {

  if (!s.solve(vars.vars)) {
    return;
  }
  std::vector<Lit> unitsLit;
  std::vector<Var> freeVars;
  s.whichAreUnits(vars.vars, unitsLit); // collect unit literals
  specs.preUpdate(unitsLit);
  for (auto i : unitsLit) {
    units_cover[i.var()] += 1.0 / (depth + 1);
  }
  std::vector<ProjVars> varConnected;
  int nbComponent =
      specs.computeConnectedComponent(varConnected, vars.vars, freeVars);
  varConnected.erase(
      std::partition(varConnected.begin(), varConnected.end(),
                     [&](const ProjVars &comp) { return comp.nbProj > 0; }),
      varConnected.end());

  std::sort(varConnected.begin(), varConnected.end(),
            [](ProjVars &a, ProjVars &b) {
              bool clean_a = a.vars.size() == a.nbProj;
              bool clean_b = b.vars.size() == b.nbProj;
              return clean_a > clean_b;
            });

  nbComponent = varConnected.size();
  if (nbComponent) {
    for (int cp = 0; cp < nbComponent; cp++) {
      ProjVars &connected = varConnected[cp];
      bool hasPriority = false;
      for (auto v : connected.vars) {
        if (specs.varIsAssigned(v) || !specs.isSelected(v))
          continue;
        if ((hasPriority = cut[v]))
          break;
      }
      if (!hasPriority) {
        continue;
      }
      Var v = score->selectVariable(connected.vars, specs);
      Lit l = Lit::makeLit(v, true);

      s.pushAssumption(l);
      probe(s, specs, score, connected, units_cover, cut, depth + 1);
      s.popAssumption();

      if (s.isInAssumption(l))
        continue;
      else if (s.isInAssumption(~l))
        probe(s, specs, score, connected, units_cover, cut, depth + 1);
      else {
        s.pushAssumption(~l);
        probe(s, specs, score, connected, units_cover, cut, depth + 1);
        s.popAssumption();
      }
    }
    specs.postUpdate(unitsLit);
    return;
  } // else we have a tautology
    //
  specs.postUpdate(unitsLit);
  return;
}

void PartitioningHeuristicStaticSingleProj::look_ahead(
    std::vector<Var> &component, std::vector<Var> &equivClass,
    std::vector<std::vector<Var>> &equivVar,
    std::vector<unsigned> &bucketNumber) {
  std::vector<int> partition(m_maxNbNodes, 0);
  std::vector<Var> considered;
  std::vector<unsigned> var2class(m_nbVar + 1, -1);

  auto get_cut = [&]() {
    std::vector<unsigned> cut;
    for (auto &e : m_hypergraph) {
      int p = partition[e.getId()];
      bool in_cut = false;
      for (auto v : e) {
        if (p != partition[v]) {
          in_cut = true;
        }
      }
      if (in_cut) {
        cut.push_back(e.getId());
      }
    }
    for (int i = cut.size() - 1; i >= 0; i--) {
      for (auto e : equivVar[var2class[cut[i]]]) {
        cut.push_back(e);
      }
    }
    for (int i = cut.size() - 1; i >= 0; i--) {
      if (m_om.isSelected(cut[i])) {
        std::swap(cut[i], cut.back());
        cut.pop_back();
      }
    }
    return cut;
  };
  struct EquivClass {
    std::vector<Var> v;
    float gain;
  };
  while (true) {
    m_hypergraphExtractor->constructHyperGraph(m_om, component, equivClass,
                                               equivVar, m_reduceFormula,
                                               considered, m_hypergraph);
    m_pm->computePartition(m_hypergraph, PartitionerManager::Level::QUALITY,
                           partition, m_multi_lvl_cuts);
    for (unsigned i = 0; i < equivVar.size(); i++) {
      for (auto v : equivVar[i]) {
        var2class[v] = i;
      }
    }
    auto cut = get_cut();

    int depth = 5;
    if (cut.size() < depth) {
      depth = cut.size();
    }
    std::vector<bool> contain(cut.size());
    std::vector<bool> hasPriority(m_nbVar + 1, false);
    std::vector<float> coverage(m_nbVar + 1);
    std::vector<EquivClass> classes;

    contain[0] = 1;
    for (int i = 1; i < depth; i++) {
      contain[i] = 1;
      // All permutation
      do {
        std::fill(hasPriority.begin(), hasPriority.end(), false);
        std::fill(coverage.begin(), coverage.end(), 0.0);
        for (int j = 0; j < cut.size(); j++) {
          if (contain[j]) {
            hasPriority[cut[j]] = true;
          }
        }
      } while (prev_permutation(contain.begin(), contain.end()));
    }
  }
}

int PartitioningHeuristicStaticSingleProj::distributePartition(
    std::vector<std::vector<unsigned>> &hypergraph, std::vector<int> &partition,
    std::vector<unsigned> &mappingEdge, std::vector<Var> &mappingVar,
    std::vector<Strata> &stack, unsigned &level, std::vector<bool> &is_proj) {
  std::vector<unsigned> cutSet, indicesFirst, indicesSecond;
  hyper_util::splitWrtPartition(m_hypergraph, partition, mappingEdge, cutSet,
                                indicesFirst, indicesSecond);

  int projV1 = 0;
  for (auto i : indicesFirst) {
    Var v = mappingVar[i];
    projV1 += is_proj[v];
  }
  int projV2 = 0;
  for (auto i : indicesSecond) {
    Var v = mappingVar[i];
    projV2 += is_proj[v];
  }

  unsigned fatherId = stack.back().fatherId;
  unsigned currentId = (cutSet.size()) ? level : fatherId;
  stack.pop_back();

  if (cutSet.size()) {
    int cut_nproj = 0;
    for (unsigned i : cutSet) {
      int v = mappingVar[i];
      cut_nproj += is_proj[v];
    }
    float ratio = float(cut_nproj) / (projV2 + projV1);
    if (ratio > m_max_cut_ratio) {
      std::vector<unsigned> all = indicesFirst;
      all.insert(all.end(), indicesSecond.begin(), indicesSecond.end());
      all.insert(all.end(), cutSet.begin(), cutSet.end());
      assignLevel(hypergraph, currentId, all, mappingVar, level);
      std::cout << "Cut to large: " << ratio << std::endl;
      return cutSet.size();
    } else {
      setCutSetBucketLevelFromEdges(hypergraph, partition, cutSet, mappingVar,
                                    level);
      assert(fatherId < m_levelInfo.size());
      m_levelInfo[fatherId].separatorLevel = level;

      level++;
      m_levelInfo.push_back({level, (unsigned)cutSet.size()});
    }
  } else {
    // special case 1.
    if (!indicesFirst.size() && indicesSecond.size()) {
      assignLevel(hypergraph, currentId, indicesSecond, mappingVar, level);
      return 0;
    }

    // special case 2.
    if (!indicesSecond.size() && indicesFirst.size()) {
      assignLevel(hypergraph, currentId, indicesFirst, mappingVar, level);
      return 0;
    }
  }

  if (projV2 > LIMIT)
    stack.push_back({currentId, indicesSecond});
  else
    assignLevel(hypergraph, currentId, indicesSecond, mappingVar, level);
  if (projV1 > LIMIT)
    stack.push_back({currentId, indicesFirst});
  else
    assignLevel(hypergraph, currentId, indicesFirst, mappingVar, level);
  return cutSet.size();
}
void PartitioningHeuristicStaticSingleProj::find_bad_vars(
    std::vector<Var> &bad) {}
void PartitioningHeuristicStaticSingleProj::computeDecomposition(
    std::vector<Var> &component, std::vector<Var> &equivClass,
    std::vector<std::vector<Var>> &equivVar,
    std::vector<unsigned> &bucketNumber) {
  using Level = PartitionerManager::Level;
  assert(m_equivClass.size() == equivClass.size());
  for (unsigned i = 0; i < equivClass.size(); i++)
    m_equivClass[i] = equivClass[i];

  // construct the hypergraph
  std::vector<Var> considered;
  m_hypergraphExtractor->constructHyperGraph(m_om, component, equivClass,
                                             equivVar, m_reduceFormula,
                                             considered, m_hypergraph);
  /*
  std::cout << "CurrentHyper Graph:" << std::endl;
  for (auto edge : m_hypergraph) {
    std::cout << "Edge: " << edge.getId()
              << " W:" << m_hypergraph.getCost()[edge.getId()] << std::endl
              << "Vert: ";
    for (auto v : edge)
      std::cout << v << " ";
    std::cout << std::endl;
  }
  */

  // save the hyper graph.
  std::vector<std::vector<unsigned>> savedHyperGraph;
  std::vector<int> weight;
  hyper_util::saveHyperGraph(savedHyperGraph, weight, m_hypergraph);

  // preparation.
  std::vector<int> partition(m_maxNbNodes, 0);
  std::vector<Var> indexToVar(m_maxNbEdges, 0);

  // init the stack with all the edges.
  std::vector<Strata> stack;
  Strata strata = {0, std::vector<unsigned>()};
  for (unsigned i = 0; i < savedHyperGraph.size(); i++)
    strata.part.push_back(i);
  stack.push_back(strata);

  // reinit the bucket for all.
  m_levelInfo.clear();
  m_levelInfo.push_back({0, 0});
  for (auto &b : m_bucketNumber)
    b = 0;
  unsigned level = 1;

  // iteratively consider sub-graph.
  std::vector<bool> is_proj(m_om.getNbVariable() + 1);
  for (Var v : component) {
    is_proj[v] = m_om.isSelected(v);
  }
  for (auto vec : equivVar) {
    bool sel = 0;
    for (auto v : vec) {
      sel |= m_om.isSelected(v);
    }
    is_proj[vec.back()] = sel;
  }
  float imbalance = 0.1;
  float ratio = float(m_om.nbSelected()) / m_om.getNbVariable();
  imbalance = std::min(imbalance + (1 - ratio), 0.8f);
  imbalance = 0.5;
  std::cout<<"imbalance set to "<<imbalance<<std::endl;
  while (stack.size()) {
    Strata &strata = stack.back();
    std::vector<unsigned> &current = strata.part;
    hyper_util::setHyperGraph(savedHyperGraph, weight, current, m_hypergraph);
    if ((m_multi_lvl_cuts > 2) && level == 1) {
      m_pm->computePartition(m_hypergraph, Level::QUALITY, partition,
                             m_multi_lvl_cuts);
      handle_multicut(m_hypergraph, savedHyperGraph, partition, current,
                      considered, stack, level);
    } else {

      m_pm->computePartition(m_hypergraph, Level::QUALITY, partition, 2,
                             imbalance);

      // get the cut and split the current set of variables.
      int cut_size = distributePartition(savedHyperGraph, partition, current,
                                         considered, stack, level, is_proj);
    }
  }

#if LOG
  if (m_log.is_open()) {
    m_log.close();
  }
#endif
  // set the equivalence.
  for (auto v : component) {
    if (m_bucketNumber[v])
      continue;

    if (v == equivClass[v])
      m_bucketNumber[v] = level;
    else
      m_bucketNumber[v] = m_bucketNumber[equivClass[v]];
  }
} // computeDecomposition

} // namespace d4
