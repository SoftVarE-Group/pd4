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
#include "utils/System.h"
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <span>
#include <unordered_map>
#include <unordered_set>

namespace d4 {

EquivExtractorProj::EquivExtractorProj(SpecManagerCnf *specs)
    : EquivExtractor(specs->getNbVariable()), m_specs(specs) {}

/*
struct Hit {
  Var v = 0;
  int cnt = 0;
};
struct SatLvl {
  std::vector<Hit> hits;
  int sat = 0;
  void hit(Var v) {
    for (auto &h : hits) {
      if (v == h.v) {
        h.cnt++;
        return;
      }
    }
    hits.push_back({v, 1});
  }
};

void count_implied_units(WrapperSolver &s, SpecManagerCnf &specs,
                         std::unordered_set<Var> &component,
                         std::vector<Var> &selection, int depth,
                         std::vector<SatLvl> &lvls) {
  std::vector<Lit> units;
  std::vector<Var> free_vars;
  static std::vector<Var> component_vec = {};
  component_vec.clear();
  component_vec.insert(component_vec.end(), component.begin(), component.end());
  if (!s.solve(component_vec)) {
    return;
  }
  lvls[depth].sat++;
  s.whichAreUnits(component_vec, units);
  specs.preUpdate(units);
  specs.freeVars(component_vec, free_vars);
  for (auto ll : units) {
    lvls[depth].hit(ll.var());
    component.erase(ll.var());
  }
  for (auto v : free_vars) {
    component.erase(v);
  }

  while (depth < selection.size() && !component.contains(selection[depth])) {
    depth++;
  }
  if (depth >= selection.size()) {
    specs.postUpdate(units);
    for (auto ll : units) {
      component.insert(ll.var());
    }
    for (auto v : free_vars) {
      component.insert(v);
    }
    return;
  }

  Lit l = Lit::makeLitTrue(selection[depth]);

  s.pushAssumption(l);
  count_implied_units(s, specs, component, selection, depth + 1, lvls);
  s.popAssumption();

  s.pushAssumption(~l);
  count_implied_units(s, specs, component, selection, depth + 1, lvls);
  s.popAssumption();

  specs.postUpdate(units);
  for (auto ll : units) {
    component.insert(ll.var());
  }
  for (auto v : free_vars) {
    component.insert(v);
  }
}

struct EqClass {
  std::vector<Var> vars;
  std::vector<Var> sel;
  int nb_proj = 0;
  int nb_nproj = 0;
  EqClass(SpecManagerCnf &specs, std::vector<Var> sel, std::vector<Var> &vars)
      : vars(vars), sel(sel) {
    for (auto v : vars) {
      nb_proj += specs.isSelected(v);
    }
    nb_nproj = vars.size() - nb_proj;
  }
  double gain() const { return (nb_nproj) / double(sel.size()); }
};

void find_implied_sets(WrapperSolver &s, SpecManagerCnf &specs,
                       std::vector<Var> &components,
                       std::vector<Var> &selection, double min_hit_rate,
                       std::vector<EqClass> &output) {
  std::vector<SatLvl> lvls(selection.size() + 1);
  std::unordered_set<Var> component_set(components.begin(), components.end());
  count_implied_units(s, specs, component_set, selection, 0, lvls);
  assert(std::unordered_set<Var>(components.begin(), components.end()) ==
         component_set);
  std::vector<Var> implied_vars;
  for (int i = 1; i <= selection.size(); i++) {
    auto &l = lvls[i];
    for (auto h : l.hits) {
      if (h.cnt / double(l.sat) >= min_hit_rate) {
        implied_vars.push_back(h.v);
      }
    }
    output.push_back(EqClass(
        specs, std::vector<Var>(selection.begin(), selection.begin() + i),
        implied_vars));
  }
}
void find_implied_sets2(WrapperSolver &s, SpecManagerCnf &specs,
                        std::vector<Var> &components,
                        std::vector<Var> &selection, double min_hit_rate,
                        std::vector<EqClass> &output) {
  for (int k = 0; k < (1 << components.size()); k++) {
  }
}
*/

void EquivExtractorProj::searchEquiv(WrapperSolver &s, std::vector<Var> &vars,
                                     std::vector<std::vector<Var>> &equivVar) {
  /*
    PrimalGraph G(*m_specs);
    int max_depth = 3;
    int max_samples = 30000;
    double max_time = 10.0;
    auto start_time = Glucose::cpuTime();
    struct VectorHash {
      size_t operator()(const std::vector<int> &v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int k = 0; k < v.size(); k++) {
          seed ^= hasher(v[k]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
      }
    };
    std::vector<EqClass> sets;
    std::vector<Var> proj_vars;
    std::vector<int> proj_freq;
    for (auto v : vars) {
      if (m_specs->isSelected(v)) {
        proj_vars.push_back(v);
        proj_freq.push_back(G.freq[v]);
        // G.freq[v1]>G.freq[v2];});
      }
    }

    for (auto &v : vars)
      m_flagVar[v] = true;
    std::vector<int> eq_prob(m_specs->getNbVariable() + 1);
    for (auto &v : vars) {
      assert((unsigned)v < m_markedVar.size());
      if (m_markedVar[v] || s.varIsAssigned(v) )
        continue;
      std::vector<Var> eqv;
      if (interCollectUnit(s, v, eqv, m_flagVar)) {
        assert(eqv.size() > 0);
        if (eqv.size() == 1)
          continue;
        sets.push_back(EqClass(*m_specs, {v}, eqv));
        for (auto &vv : eqv) {
          m_markedVar[vv] = true;
        }
      }
    }
    if (proj_vars.size() <= 10) {
      return;
    }
    std::unordered_set<std::vector<Var>, VectorHash> seen;
    std::mt19937 gen(42);
    // std::discrete_distribution<int> dist(proj_freq.begin(),proj_freq.end());
    std::uniform_int_distribution<int> dist(0, proj_vars.size() - 1);

    while ((Glucose::cpuTime() - start_time) < max_time &&
           sets.size() < max_samples && false) {
      Var start = proj_vars[dist(gen)];
      std::vector<Var> sel = {start};
      auto &edges = G.edges[start];
      while (sel.size() < max_depth) {
        std::uniform_int_distribution<int> edist(0, edges.size() - 1);
        Var next = edges[edist(gen)];
        if (std::find(sel.begin(), sel.end(), next) != sel.end()) {
          break;
        } else {
          edges = G.edges[next];
          sel.push_back(next);
        }
      }
      while (sel.size() < max_depth) {
        Var next = dist(gen);
        if (std::find(sel.begin(), sel.end(), next) == sel.end()) {
          edges = G.edges[next];
          sel.push_back(next);
        }
      }
      if (seen.insert(sel).second) {
        find_implied_sets(s, *m_specs, vars, sel, 0.99, sets);
      }
    }

    sort(sets.begin(), sets.end(), [&](const EqClass &lhs, const EqClass &rhs) {
      return lhs.gain() > rhs.gain();
    });
    std::vector<bool> cover(m_specs->getNbVariable() + 1, false);
    for (auto &i : sets) {
      bool ok = true;
      for (auto v : i.vars) {
        if (cover[v]) {
          ok = false;
          break;
        }
      }
      if (i.vars.size() <= i.sel.size() || i.nb_nproj == 0) {
        continue;
      }
      if (ok) {
        for (auto v : i.vars) {
          cover[v] = true;
        }
        std::cout << "Classs: "
                  << "[";
        for (auto v : i.sel) {
          std::cout << v << ",";
        }
        std::cout << "]:[";
        for (auto v : i.vars) {
          std::cout << v << ",";
        }
        std::cout << "]" << std::endl;
        equivVar.push_back(i.vars);
      }
    }
  */

  std::vector<Var> reinit;
  for (auto &v : vars)
    m_flagVar[v] = true;
  for (auto &v : vars) {
    assert((unsigned)v < m_markedVar.size());
    if (m_markedVar[v] || s.varIsAssigned(v) || !m_specs->isSelected(v))
      continue;

    std::vector<Var> eqv;
    if (interCollectUnit(s, v, eqv, m_flagVar)) {
      assert(eqv.size() > 0);
      
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

  for (int i = equivVar.size() - 1; i >= 0; i--) {
    std::cout << "Classs: "
              << "[";
    for (auto v : equivVar[i]) {
      std::cout << v << ",";
    }
    std::cout << "]:[";
    for (auto v : equivVar[i]) {
      std::cout << v << ",";
    }
  }

} // searchEquiv

} // namespace d4
