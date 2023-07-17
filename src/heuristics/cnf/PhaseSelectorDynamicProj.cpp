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
#include "PhaseSelectorDynamicProj.hpp"

#include <ostream>

namespace d4 {

/**
   Constructor.

   @param[in] limitPhase, give the limit number of variables before switching.
*/
PhaseSelectorDynamicProj::PhaseSelectorDynamicProj(
    PartitioningHeuristicStaticSingle *staticPartitioner, double limitRatio,
    std::ostream &out)
    : PhaseSelectorManager(staticPartitioner) {
  out << "c [CONSTRUCTOR] Switching between static and dynamic decomposition:"
      << " dynamic proj(" << limitRatio << ")\n";

  m_limitRatio = limitRatio;
} // constructor

/**
   Check out if the current tree decomposition is still OK. To do it we check
   out if the current decomposition is not too unbalanced.

   @param[in] component, the current set of variables.
 */
bool PhaseSelectorDynamicProj::isStillOk(std::vector<Var> &component) {
  if (!component.size() || component.size() < 100)
    return true;
  std::vector<Var> cut;
  m_staticPartitioner->computeCutSet(component, cut);
  for (auto v : cut) {
    if (m_staticPartitioner->specs().isSelected(v)) {
      return true;
    }
  }
  return false;
} // isStillok

} // namespace d4
