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

#include "EquivExtractorProj.hpp"
#include "UnionFind.hpp"
#include "src/utils/UnionFind.hpp"
#include <boost/dynamic_bitset.hpp>
#include <chrono>
#include <random>
#include <span>
#include <unordered_map>
#include <unordered_set>

namespace d4 {
EquivExtractorProj::EquivExtractorProj(SpecManagerCnf *specs)
    : EquivExtractor(specs->getNbVariable()), m_specs(specs) {}



void EquivExtractorProj::searchEquiv(WrapperSolver &s, std::vector<Var> &vars,
                                     std::vector<std::vector<Var>> &equivVar) {
  std::vector<Var> reinit;
  for (auto &v : vars)
    m_flagVar[v] = true;
  for (auto &v : vars) {
    assert((unsigned)v < m_markedVar.size());
    if (m_markedVar[v] || s.varIsAssigned(v)||m_specs->isSelected(v))
      continue;

    std::vector<Var> eqv;
    if (interCollectUnit(s, v, eqv, m_flagVar)) {
      assert(eqv.size() > 0);
      if (eqv.size() == 1)
        continue;
      equivVar.push_back(eqv);
      for (auto &vv : eqv) {
        m_markedVar[vv] = true;
        reinit.push_back(vv);
      }
    }
  }

  for (auto &v : reinit)
    m_markedVar[v] = false;
  for (auto &v : vars)
    m_flagVar[v] = false;

  // fusion the equivalence classes that share variables.
  for (unsigned i = 0; i < equivVar.size(); i++) {
    for (auto &v : equivVar[i])
      m_markedVar[v] = true;

    unsigned j = i + 1;
    while (j < equivVar.size()) {
      bool share = false;
      for (auto &v : equivVar[j]) {
        share = m_markedVar[v];
        if (share)
          break;
      }

      if (!share)
        j++;
      else {
        for (auto &v : equivVar[j]) {
          if (!m_markedVar[v]) {
            m_markedVar[v] = true;
            equivVar[i].push_back(v);
          }
        }

        equivVar[j].clear();
        j = i + 1;
      }
    }

    for (auto &v : equivVar[i])
      m_markedVar[v] = false;
  }

  // remove the empty list
  unsigned j = 0;
  for (unsigned i = 0; i < equivVar.size(); i++) {
    if (equivVar[i].size()) {
      if (i != j)
        equivVar[j] = equivVar[i];
      j++;
    }
  }
  equivVar.resize(j);

} // searchEquiv

} // namespace d4
