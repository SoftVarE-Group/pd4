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

#include "ScoringMethodVsads2.hpp"

namespace d4 {

/**
   Constructor.

   @param[in] o, the specification of a CNF problem.
   @param[in] a, an activity manager linked to a solver.
 */
ScoringMethodVsads2::ScoringMethodVsads2(SpecManagerCnf &o, ActivityManager &a)
    : om(o), activity(a) {} // constructor

/**
   The well-known VSADS heuristic.

   Sang, T.; Beame, P.; and Kautz, H. Heuristics for Fast
   Exact Model Counting. In Proceedings of the 8th International
   Conference on Theory and Applications of Satisfiability Testing.

   @param[in] v, the variable we want the score.
 */
double ScoringMethodVsads2::computeScore(Var v) {
  throw std::runtime_error("Unsupported operation");
  return 0;
}

Var ScoringMethodVsads2::selectVariable(std::vector<Var> &vars,
                                        std::function<bool(Var)> can_select) {

  Var ret = var_Undef;
  double bestScore_act = -1;
  int bestScore_freq = -1;
  for (auto &v : vars) {
    if (!can_select(v))
      continue;
    int freq = om.getNbClause(v);
    if (ret == var_Undef || freq > bestScore_freq) {
      ret = v;
      bestScore_freq = freq;
      bestScore_act = activity.getActivity(v);
    } else if (freq == bestScore_freq) {
      int act = activity.getActivity(v);
      if (act > bestScore_act) {
        ret = v;
        bestScore_act = act;
      }
    }
  }
  return ret;
}

void ScoringMethodVsads2::decay(){
    activity.decay();

}

} // namespace d4
