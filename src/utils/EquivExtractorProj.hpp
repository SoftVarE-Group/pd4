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
#include "src/utils/UnionFind.hpp"

namespace d4 {
class EquivExtractorProj : public EquivExtractor {
private:
  SpecManagerCnf *m_specs = 0;
  UnionFind m_uf;
  std::vector<unsigned int> m_idx_clauses;
  std::vector<std::vector<Var>> m_bags;
  std::vector<int> m_active_bags;

  bool interCollectUnit(WrapperSolver &s, Var v, std::vector<Var> &listVarPU,
                        std::vector<bool> &flagVar);

public:
  EquivExtractorProj() = default; // empty constructor
  EquivExtractorProj(SpecManagerCnf *specs);

  void searchEquiv(WrapperSolver &s, std::vector<Var> &v,
                   std::vector<std::vector<Var>> &equivVar) final;
};
} // namespace d4
