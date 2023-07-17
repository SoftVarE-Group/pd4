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
#include "PartitionerPatohProj.hpp"

#include <iostream>
#include <vector>

#include "3rdParty/patoh/patoh.h"
#include "src/exceptions/OptionException.hpp"

namespace d4 {

/**
   Constructor.

   @param[in] maxNodes, the maximal number of nodes.
   @param[in] maxEdges, the maximal number of hyper edges.
 */
PartitionerPatohProj::PartitionerPatohProj(unsigned maxNodes, unsigned maxEdges,
                                           unsigned maxSumEdgeSize,
                                           std::ostream &out)
    : PartitionerPatoh(maxNodes, maxEdges, maxSumEdgeSize, out) {
  m_cost = new int[maxEdges];

} // constructor

/**
   Destructor.
 */
PartitionerPatohProj::~PartitionerPatohProj() { delete[] m_cost; } // destructor

/**
   Get a partition from the hypergraph.

   @param[in] hypergraph, the graph we search for a partition.
   @param[out] parition, the resulting partition (we suppose it is allocated).
 */
void PartitionerPatohProj::computePartition(HyperGraph &hypergraph, Level level,
                                            std::vector<int> &partition) {
  std::vector<unsigned> elts;

  // graph initialization and shift the hypergraph
  unsigned sizeXpins = 0;
  int posPins = 0;

  for (auto &edge : hypergraph) {
    m_cost[sizeXpins] = hypergraph.getCost()[edge.getId()];
    m_xpins[sizeXpins++] = posPins;
    for (auto x : edge) {
      assert(x < m_markedNodes.size());
      if (!m_markedNodes[x]) {
        m_markedNodes[x] = true;
        m_mapNodes[x] = elts.size();
        elts.push_back(x);
      }

      m_pins[posPins++] = m_mapNodes[x];
    }
  }

  for (auto &x : elts)
    m_markedNodes[x] = false;
  m_xpins[sizeXpins] = posPins;

  // hypergraph partitioner
  PaToH_Parameters args;
  switch (level) {
  case NORMAL:
    PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
    break;
  case SPEED:
    PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_SPEED);
    break;
  case QUALITY:
    PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_QUALITY);
    break;
  default:
    throw(OptionException("Wrong option given to the partioner.", __FILE__,
                          __LINE__));
  }

  args._k = 2;
  args.seed = 2911;

  int cut;
  PaToH_Alloc(&args, elts.size(), sizeXpins, 1, m_cwghts, m_cost, m_xpins,
              m_pins);
  PaToH_Part(&args, elts.size(), sizeXpins, 1, 0, m_cwghts, m_cost, m_xpins,
             m_pins, NULL, m_partvec, m_partweights, &cut);

  for (unsigned i = 0; i < elts.size(); i++)
    partition[elts[i]] = m_partvec[i];
  PaToH_Free();
} // computePartition

} // namespace d4
