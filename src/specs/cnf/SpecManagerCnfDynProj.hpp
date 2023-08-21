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
#include <cassert>
#include <src/problem/ProblemManager.hpp>
#include <src/problem/ProblemTypes.hpp>
#include <vector>

#include "SpecManagerCnfDyn.hpp"

struct ProjClauseInfo{
    int nbProj;
    int nbNProj;
};
namespace d4 {
class SpecManagerCnfDynProj : public SpecManagerCnfDyn {
 private:
    std::vector<bool> m_is_mixed;
    std::vector<ProjClauseInfo> m_clause_proj_vars;
 public:
  SpecManagerCnfDynProj(ProblemManager &p);
  virtual ~SpecManagerCnfDynProj();
  void preUpdate(std::vector<Lit> &lits);
  void postUpdate(std::vector<Lit> &lits);

  // we cannot use this function here
  inline void initialize(std::vector<Var> &setOfVar, std::vector<Lit> &units) {
    assert(0);
  }
};
}  // namespace d4
