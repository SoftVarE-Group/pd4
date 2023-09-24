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
 * You should have received a copy of the GNU Lesser General Public LicenseProj
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <ostream>

#include "PartitioningHeuristicStaticSingle.hpp"
#include "PhaseSelectorManager.hpp"
#include "src/heuristics/ScoringMethod.hpp"
#include "src/solvers/WrapperSolver.hpp"

namespace d4 {
class PhaseSelectorManager;
namespace hyper_util {
void saveHyperGraph(std::vector<std::vector<unsigned>> &savedHyperGraph,
                    std::vector<int> &savedCost, HyperGraph &h);

void setHyperGraph(std::vector<std::vector<unsigned>> &savedHyperGraph,
                   std::vector<int> &savedCost, std::vector<unsigned> &indices,
                   HyperGraph &hypergraph);

} // namespace hyper_util

class PartitioningHeuristicStaticSingleProj
    : public PartitioningHeuristicStaticSingle {
protected:
  int m_multi_lvl_cuts = 2;
  ScoringMethod *m_score;

  void distributePartition(std::vector<std::vector<unsigned>> &hypergraph,
                           std::vector<int> &partition,
                           std::vector<unsigned> &mappingEdge,
                           std::vector<Var> &mappingVar,
                           std::vector<Strata> &stack, unsigned &level);

  void handle_multicut(HyperGraph &hypergraph,
                       std::vector<std::vector<unsigned>> &hypergraph_list,
                       std::vector<int> &partition,
                       std::vector<unsigned> &mappingEdge,
                       std::vector<Var> &mappingVar, std::vector<Strata> &stack,
                       unsigned &level);

public:
  PartitioningHeuristicStaticSingleProj(po::variables_map &vm, WrapperSolver &s,
                                        SpecManager &om, std::ostream &out);

  PartitioningHeuristicStaticSingleProj(po::variables_map &vm, WrapperSolver &s,
                                        SpecManager &om, int nbClause,
                                        int nbVar, int sumSize,
                                        std::ostream &out);
  void handle_multicut();
  void look_ahead(std::vector<Var> &component, std::vector<Var> &equivClass,
                  std::vector<std::vector<Var>> &equivVar,
                  std::vector<unsigned> &bucketNumber);
  void find_bad_vars(std::vector<Var> &bad);
  void computeDecomposition(std::vector<Var> &component,
                            std::vector<Var> &equivClass,
                            std::vector<std::vector<Var>> &equivVar,
                            std::vector<unsigned> &bucketNumber) override;
};
} // namespace d4
