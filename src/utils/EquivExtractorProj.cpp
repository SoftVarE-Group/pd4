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

namespace d4 {
EquivExtractorProj::EquivExtractorProj(SpecManagerCnf *specs)
    : EquivExtractor(specs->getNbClause()), m_specs(specs),
      m_uf(specs->getNbVariable() + 1), m_bags(specs->getNbVariable() + 1) {}

bool EquivExtractorProj::interCollectUnit(WrapperSolver &s, Var v,
                                          std::vector<Var> &listVarPU,
                                          std::vector<bool> &flagVar) {
  std::vector<Lit> listVarPosLit, listVarNegLit;
  if (!s.decideAndComputeUnit(Lit::makeLit(v, false), listVarPosLit))
    return false;
  if (!s.decideAndComputeUnit(Lit::makeLit(v, true), listVarNegLit))
    return false;

  // intersection.
  for (auto &l : listVarPosLit)
    if (flagVar[l.var()])
      m_markedVarInter[l.var()] = true;
  for (auto &l : listVarNegLit)
    if (m_markedVarInter[l.var()] && m_specs->isSelected(v))
      listVarPU.push_back(l.var());
  for (auto &l : listVarPosLit)
    m_markedVarInter[l.var()] = false;

  return true;
} // interCollectUnit

/**
   Research equivalences in the set of variable v.

   @param[in] s, a wrapper to a solver.
   @param[in] v, the set of variables we search in.
   @param[out] equivVar, le resulting equivalences.
 */
void EquivExtractorProj::searchEquiv(WrapperSolver &s, std::vector<Var> &vars,
                                     std::vector<std::vector<Var>> &equivVar) {
  std::vector<Var> reinit;

  for (auto &v : vars)
    m_flagVar[v] = true;
  {
    m_specs->getCurrentClauses(m_idx_clauses, vars);
    for (auto i : m_idx_clauses) {
      Var root = var_Undef;
      for (auto l : m_specs->getClause(i)) {
        if (m_specs->isSelected(l.var()) || s.varIsAssigned(l.var()))
          continue;
        if (root != var_Undef) {
          m_uf.union_sets(root, l.var());
        }
        root = l.var();
      }
    }
    for (auto v : vars) {
      if (m_specs->isSelected(v) || s.varIsAssigned(v))
        continue;
      reinit.push_back(v);
      Var root = m_uf.find_set(v);
      if (m_uf.size[root] > 1) {
        m_bags[root].push_back(v);
        if (!m_markedVar[root]) {
          m_active_bags.push_back(root);
          m_markedVar[root] = true;
          reinit.push_back(root);
        }
      }
    }
    for (auto i : m_active_bags) {
      equivVar.push_back(std::move(m_bags[i]));
    }
    m_uf.clear();
  }
  for (auto &v : vars) {
    assert((unsigned)v < m_markedVar.size());
    if (m_markedVar[v] || s.varIsAssigned(v) || !m_specs->isSelected(v))
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
