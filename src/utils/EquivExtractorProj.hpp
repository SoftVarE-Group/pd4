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
#pragma once
#include "EquivExtractor.hpp"
#include "src/hyperGraph/HyperGraphExtractor.hpp"
#include "src/partitioner/PartitionerManager.hpp"
#include "src/utils/UnionFind.hpp"
#include <random>

namespace d4 {
struct PrimalGraph {
  std::vector<std::vector<int>> edges;
  std::vector<std::vector<bool>> connected;
  std::vector<int> freq, bin_freq;
  PrimalGraph(SpecManagerCnf &specs)
      : edges(specs.getNbVariable() + 1),
        connected(specs.getNbVariable() + 1,
                  std::vector<bool>(specs.getNbVariable() + 1, false)),
        freq(specs.getNbVariable() + 1), bin_freq(specs.getNbVariable() + 1) {

    for (int k = 0; k < specs.getNbClause(); k++) {
      auto &clause = specs.getClause(k);

      for (int i = 0; i < clause.size(); i++) {
        Var v1 = clause[i].var();
        if (!specs.isSelected(v1))
          continue;
        freq[v1]++;
        if (clause.size() == 2) {
          bin_freq[v1]++;
        }
        for (int j = i + 1; j < clause.size(); j++) {
          Var v2 = clause[j].var();

          if (!specs.isSelected(v2))
            continue;
          if (connected[v1][v2])
            continue;
          edges[v1].push_back(v2);
          edges[v2].push_back(v1);
          connected[v1][v2] = true;
          connected[v2][v1] = true;
        }
      }
    }
  }
};
class EquivExtractorProj : public EquivExtractor {
private:
  SpecManagerCnf *m_specs = 0;

public:
  EquivExtractorProj() = default; // empty constructor
  EquivExtractorProj(SpecManagerCnf *specs);

  void count_implied_units(WrapperSolver &s, const std::vector<Var> &vars,
                           std::vector<Var> &units, double branches);
  void searchEquiv(WrapperSolver &s, std::vector<Var> &v,
                   std::vector<std::vector<Var>> &equivVar) final;
};
} // namespace d4
