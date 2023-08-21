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

#include "SpecManagerCnfDynProj.hpp"

namespace d4 {

/**
   OccurrenceManager constructor. This function initialized the
   structures used.

   @param[in] _clauses, the set of clauses
   @param[in] _nbVar, the number of variables in the problem
 */
SpecManagerCnfDynProj::SpecManagerCnfDynProj(ProblemManager &p)
    : SpecManagerCnfDyn(p) {
  m_clause_proj_vars.resize(getNbClause());
  m_is_mixed.resize(getNbClause(), false);
  for (int i = 0; i < m_clauses.size(); i++) {
    auto &cl = m_clauses[i];
    int nb_proj = 0;
    int nb_nproj = 0;
    for (auto v : cl) {
      if (isSelected(v.var())) {
        nb_proj++;
      } else {
        nb_nproj++;
      }
    }
    if (nb_nproj != 0 && nb_proj != 0) {
      m_is_mixed[i] = true;
      for (auto v : cl) {
        m_occurrence[v.intern()].nbMixed++;
      }
    }
    m_clause_proj_vars[i] = {nb_proj, nb_nproj};
  }

} // SpecManagerCnfDyn
  //
SpecManagerCnfDynProj::~SpecManagerCnfDynProj() {
#ifdef DEBUG
  std::vector<int> mixed_check(m_occurrence.size());
  for (int i = 0; i < m_clauses.size(); i++) {
    auto &cl = m_clauses[i];
    int nb_proj = 0;
    int nb_nproj = 0;
    for (auto v : cl) {
      if (isSelected(v.var())) {
        nb_proj++;
      } else {
        nb_nproj++;
      }
    }
    if (nb_nproj != 0 && nb_proj != 0) {
      m_is_mixed[i] = true;
      for (auto v : cl) {
        mixed_check[v.intern()]++;
      }
    }
  }
  for (int i = 0; i < mixed_check.size(); i++) {
    assert(mixed_check[i] == m_occurrence[i].nbMixed);
  }
  for (auto &i : m_infoClauses) {
    assert(i.nbUnsatProj == 0);
    assert(i.nbUnsat == 0);
    assert(i.nbSat == 0);
  }

#endif
}

/**
   Update the occurrence list w.r.t. a new set of assigned variables.
   It's important that the order is conserved between the moment where
   we assign and tVsadshe moment we unassign.

   @param[in] lits, the new assigned variables
 */
void SpecManagerCnfDynProj::preUpdate(std::vector<Lit> &lits) {
  m_reviewWatcher.resize(0);
  auto removed_mixed_proj = [&](int cl) {
    int unsat_proj = m_infoClauses[cl].nbUnsatProj;
    int unsat_nproj = m_infoClauses[cl].nbUnsat - unsat_proj;
    return (unsat_proj == m_clause_proj_vars[cl].nbProj ||
            unsat_nproj != m_clause_proj_vars[cl].nbNProj) &&
           m_is_mixed[cl];
  };
  auto removed_mixed_nproj = [&](int cl) {
    int unsat_proj = m_infoClauses[cl].nbUnsatProj;
    int unsat_nproj = m_infoClauses[cl].nbUnsat - unsat_proj;
    return (unsat_proj != m_clause_proj_vars[cl].nbProj ||
            unsat_nproj == m_clause_proj_vars[cl].nbNProj) &&
           m_is_mixed[cl];
  };
  auto removed_mixed = [&](int cl, Var v) {
    if (isSelected(v)) {
      return removed_mixed_proj(cl);
    } else {
      return removed_mixed_nproj(cl);
    }
  };
  for (auto &l : lits) {
    m_currentValue[l.var()] = l.sign() ? l_False : l_True;
    // not binary clauses.
    for (unsigned i = 0; i < m_occurrence[l.intern()].nbNotBin; i++) {
      int idxCl = m_occurrence[l.intern()].notBin[i];
      m_infoClauses[idxCl].nbSat++;
      for (auto &ll : m_clauses[idxCl])
        if (m_currentValue[ll.var()] == l_Undef) {
          m_occurrence[ll.intern()].removeNotBin(idxCl);
          m_occurrence[ll.intern()].nbMixed -= m_is_mixed[idxCl];
          assert(m_occurrence[ll.intern()].nbMixed >= 0);
        }
    }

    for (unsigned i = 0; i < m_occurrence[(~l).intern()].nbNotBin; i++) {
      int idxCl = m_occurrence[(~l).intern()].notBin[i];
      m_infoClauses[idxCl].nbUnsat++;
      m_infoClauses[idxCl].nbUnsatProj += isSelected(l.var());
      if (m_infoClauses[idxCl].watcher == ~l)
        m_reviewWatcher.push_back(idxCl);
      if (removed_mixed(idxCl, l.var())) {
        for (auto ll : m_clauses[idxCl]) {
          if (!isSelected(ll.var())) {
            break;
          }
          m_occurrence[ll.intern()].nbMixed -= 1;
          assert(m_occurrence[ll.intern()].nbMixed >= 0);
        }
      }
    }

    // binary clauses.
    for (unsigned i = 0; i < m_occurrence[l.intern()].nbBin; i++) {
      int idxCl = m_occurrence[l.intern()].bin[i];
      m_infoClauses[idxCl].nbSat++;
      for (auto &ll : m_clauses[idxCl])
        if (m_currentValue[ll.var()] == l_Undef) {
          m_occurrence[ll.intern()].removeBin(idxCl);
          m_occurrence[ll.intern()].nbMixed -= m_is_mixed[idxCl];
          assert(m_occurrence[ll.intern()].nbMixed >= 0);
        }
    }
    for (unsigned i = 0; i < m_occurrence[(~l).intern()].nbBin; i++) {
      int idxCl = m_occurrence[(~l).intern()].bin[i];
      m_infoClauses[idxCl].nbUnsat++;
      m_infoClauses[idxCl].nbUnsatProj += isSelected(l.var());
      if (m_infoClauses[idxCl].watcher == ~l)
        m_reviewWatcher.push_back(idxCl);
      if (removed_mixed(idxCl, l.var())) {
        for (auto ll : m_clauses[idxCl]) {
          if (!isSelected(ll.var())) {
            break;
          }
          m_occurrence[ll.intern()].nbMixed -= 1;
          assert(m_occurrence[ll.intern()].nbMixed >= 0);
        }
      }
    }
  }

  // we search another non assigned literal if requiered
  for (auto &idxCl : m_reviewWatcher) {
    if (m_infoClauses[idxCl].nbSat)
      continue;

    for (auto &l : m_clauses[idxCl]) {
      if (m_currentValue[l.var()] == l_Undef) {
        m_infoClauses[idxCl].watcher = l;
        break;
      }
    }
  }

} // preUpdate

/**
   Update the occurrence list w.r.t. a new set of unassigned variables.
   It's important that the order is conserved between the moment where
   we assign and the moment we unassign.

   @param[in] lits, the new assigned variables
 */
void SpecManagerCnfDynProj::postUpdate(std::vector<Lit> &lits) {
  auto removed_mixed_proj = [&](int cl) {
    int unsat_proj = m_infoClauses[cl].nbUnsatProj;
    int unsat_nproj = m_infoClauses[cl].nbUnsat - unsat_proj;
    return (unsat_proj == m_clause_proj_vars[cl].nbProj ||
            unsat_nproj != m_clause_proj_vars[cl].nbNProj) &&
           m_is_mixed[cl];
  };
  auto removed_mixed_nproj = [&](int cl) {
    int unsat_proj = m_infoClauses[cl].nbUnsatProj;
    int unsat_nproj = m_infoClauses[cl].nbUnsat - unsat_proj;
    return (unsat_proj != m_clause_proj_vars[cl].nbProj ||
            unsat_nproj == m_clause_proj_vars[cl].nbNProj) &&
           m_is_mixed[cl];
  };
  auto removed_mixed = [&](int cl, Var v) {
    if (isSelected(v)) {
      return removed_mixed_proj(cl);
    } else {
      return removed_mixed_nproj(cl);
    }
  };
  for (int i = lits.size() - 1; i >= 0; i--) {
    Lit l = lits[i];
    // for the no binary clauses.
    for (unsigned i = 0; i < m_occurrence[l.intern()].nbNotBin; i++) {
      int idxCl = m_occurrence[l.intern()].notBin[i];
      m_infoClauses[idxCl].nbSat--;
      assert(!m_infoClauses[idxCl].nbSat);

      for (auto &ll : m_clauses[idxCl])
        if (m_currentValue[ll.var()] == l_Undef) {
          m_occurrence[ll.intern()].addNotBin(idxCl);
          m_occurrence[ll.intern()].nbMixed += m_is_mixed[idxCl];
        }
    }

    for (unsigned i = 0; i < m_occurrence[(~l).intern()].nbNotBin; i++) {
      int idxCl = m_occurrence[(~l).intern()].notBin[i];
      m_infoClauses[idxCl].nbUnsat--;
      m_infoClauses[idxCl].nbUnsatProj -= isSelected(l.var());
      if (removed_mixed(idxCl, l.var())) {
        for (auto ll : m_clauses[idxCl]) {
          if (!isSelected(ll.var())) {
            break;
          }
          m_occurrence[ll.intern()].nbMixed += 1;
        }
      }
    }

    // for the binary clauses.
    for (unsigned i = 0; i < m_occurrence[l.intern()].nbBin; i++) {
      int idxCl = m_occurrence[l.intern()].bin[i];
      m_infoClauses[idxCl].nbSat--;
      assert(!m_infoClauses[idxCl].nbSat);

      for (auto &ll : m_clauses[idxCl])
        if (m_currentValue[ll.var()] == l_Undef) {
          m_occurrence[ll.intern()].addBin(idxCl);
          m_occurrence[ll.intern()].nbMixed += m_is_mixed[idxCl];
        }
    }

    for (unsigned i = 0; i < m_occurrence[(~l).intern()].nbBin; i++) {
      int idxCl = m_occurrence[(~l).intern()].bin[i];
      m_infoClauses[idxCl].nbUnsat--;
      m_infoClauses[idxCl].nbUnsatProj -= isSelected(l.var());
      if (removed_mixed(idxCl, l.var())) {
        for (auto ll : m_clauses[idxCl]) {
          if (!isSelected(ll.var())) {
            break;
          }
          m_occurrence[ll.intern()].nbMixed += 1;
        }
      }
    }

    m_currentValue[l.var()] = l_Undef;
  }
} // postUpdate

} // namespace d4
