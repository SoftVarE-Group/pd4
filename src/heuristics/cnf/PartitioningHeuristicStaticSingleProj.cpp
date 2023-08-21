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

namespace d4 {
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
    

    } // constructor

/**
   Save the current hyper graph.

   @param[out] savedHyperGraph, the structure where is saved the graph.
*/
void PartitioningHeuristicStaticSingleProj::saveHyperGraph(
    std::vector<std::vector<unsigned>> &savedHyperGraph,
    std::vector<int> &savedCost) {
  for (auto edge : m_hypergraph) {
    savedHyperGraph.push_back(std::vector<unsigned>());
    std::vector<unsigned> &tmp = savedHyperGraph.back();
    for (auto v : edge)
      tmp.push_back(v);
    savedCost.push_back(m_hypergraph.getCost()[edge.getId()]);
  }
} // savedHyperGraph

/**
   Set the hyper graph regarding the given set of variables and the saved
   hyper graph.

   @param[in] savedHyperGraph, the current hyper graph.
   @param[in] indices, the current set of edges' indices.
   @param[out] hypergraph, the computed hyper graph.
*/
void PartitioningHeuristicStaticSingleProj::setHyperGraph(
    std::vector<std::vector<unsigned>> &savedHyperGraph,
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

    //hypergraph.getCost()[hypergraph.getSize()] =
    //    savedCost[idxEdge] > 1 ? std::min<int>(2, indices.size()*0.1) : 1;
    hypergraph.incSize();
  }
}
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
  saveHyperGraph(savedHyperGraph, weight);

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
  while (stack.size()) {
    Strata &strata = stack.back();
    std::vector<unsigned> &current = strata.part;
    int hasProj = 0;
    for (auto i : current) {
      Var v = considered[i];
      hasProj += is_proj[v];
      if(hasProj>3){
          break;
      }
    }
    if (hasProj<=3) {
      setBucketLevelFromEdges(savedHyperGraph, current, considered,level);
      stack.pop_back();
      continue;
    }
    setHyperGraph(savedHyperGraph, weight, current, m_hypergraph);

    m_pm->computePartition(m_hypergraph, Level::QUALITY, partition);

    // get the cut and split the current set of variables.
    distributePartition(savedHyperGraph, partition, current, considered, stack,
                        level);
  }

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
