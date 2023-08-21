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

#include "ScoringMethodPVsads.hpp"

namespace d4 {

/**
   Constructor.

   @param[in] o, the specification of a CNF problem.
   @param[in] a, an activity manager linked to a solver.
 */
ScoringMethodPVsads::ScoringMethodPVsads(SpecManagerCnf &o, ActivityManager &a,
                                         int x, int y, int z)
    : om(o), activity(a), x(x), y(y), z(z) {} // constructor

/**
   The well-known VSADS heuristic.

   Sang, T.; Beame, P.; and Kautz, H. Heuristics for Fast
   Exact Model Counting. In Proceedings of the 8th International
   Conference on Theory and Applications of Satisfiability Testing.

   @param[in] v, the variable we want the score.
 */
double ScoringMethodPVsads::computeScore(Var v) {
  //std::cout << om.getNbMixedClause(v) << std::endl;

  return activity.getCountConflict(v) * x + om.getNbClause(v) * y +
         om.getNbMixedClause(v) * z;
}

} // namespace d4
