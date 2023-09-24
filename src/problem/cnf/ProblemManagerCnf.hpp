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

#include "../ProblemManager.hpp"
#include "../ProblemTypes.hpp"
#include "src/solvers/WrapperSolver.hpp"
#include "src/utils/IDIDFunc.hpp"

namespace d4 {
class ProblemManagerCnf : public ProblemManager {
protected:
  std::vector<std::vector<Lit>> m_clauses;
  std::vector<std::vector<Lit>> m_learnt;
  IDIDFunc m_projMap;
  int m_freevars;
public:
  ProblemManagerCnf();
  ProblemManagerCnf(int nbVar, std::vector<double> &weightLit,
                    std::vector<double> &weightVar, std::vector<Var> &selected,int freevars=0);
  ProblemManagerCnf(ProblemManager *problem);
  ProblemManagerCnf(std::string &nameFile, std::string &proj_vars);
  ~ProblemManagerCnf();
  void normalize() override;
  void normalizeInner() override;

  void display(std::ostream &out) override;
  std::vector<std::vector<Lit>> &getClauses() { return m_clauses; }
  std::vector<std::vector<Lit>> &getLearnt() { return m_learnt; }
  IDIDFunc &getProjMap() { return m_projMap; }
  void setClauses(std::vector<std::vector<Lit>> &clauses) {
    m_clauses = clauses;
  }
  void displayStat(std::ostream &out, std::string startLine) override;
  ProblemManager *getUnsatProblem() override;
  ProblemManager *getConditionedFormula(std::vector<Lit> &units) override;
  int freeVars() override;
};
} // namespace d4
