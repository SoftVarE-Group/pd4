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
#include <boost/dynamic_bitset.hpp>
#include <chrono>
#include <random>
#include <span>
#include <unordered_map>
#include <unordered_set>

namespace d4 {
EquivExtractorProj::EquivExtractorProj(SpecManagerCnf *specs)
    : EquivExtractor(specs->getNbVariable()), m_specs(specs) {}

/**
   Research equivalences in the set of variable v.

   @param[in] s, a wrapper to a solver.
   @param[in] v, the set of variables we search in.
   @param[out] equivVar, le resulting equivalences.
 */
template <class D, class W, class URBG>
void weighted_shuffle(D first, D last, W first_weight, W last_weight,
                      URBG &&g) {
  while (first != last and first_weight != last_weight) {
    std::discrete_distribution dd(first_weight, last_weight);
    auto i = dd(g);
    if (i) {
      std::iter_swap(first, std::next(first, i));
      std::iter_swap(first_weight, std::next(first_weight, i));
    }
    ++first;
    ++first_weight;
  }
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  out << "{";
  size_t last = v.size() - 1;
  for (size_t i = 0; i < v.size(); ++i) {
    out << v[i];
    if (i != last)
      out << ", ";
  }
  out << "}";
  return out;
}

void search_units(WrapperSolver &s, std::vector<Var> &vars,
                  std::vector<Var> &sel, int depth,
                  std::vector<std::vector<int>> &hits,
                  std::vector<int> &sat_cnt) {
  if (!s.solve(vars)) {
    return;
  }
  sat_cnt[depth]++;
  std::vector<Lit> units;
  s.whichAreUnits(vars, units);
  for (auto u : units) {
    hits[depth][u.var()]++;
  }
  if (depth < sel.size()) {
    Var v = sel[depth];
    Lit l = Lit::makeLitTrue(v);
    if (s.isInAssumption(l)) {
      search_units(s, vars, sel, depth + 1, hits, sat_cnt);
    } else {
      s.pushAssumption(l);
      search_units(s, vars, sel, depth + 1, hits, sat_cnt);
      s.popAssumption();
    }
    if (s.isInAssumption(~l)) {
      search_units(s, vars, sel, depth + 1, hits, sat_cnt);
    } else {
      s.pushAssumption(~l);
      search_units(s, vars, sel, depth + 1, hits, sat_cnt);
      s.popAssumption();
    }
  }
}

void EquivExtractorProj::searchEquiv(WrapperSolver &s, std::vector<Var> &vars,
                                     std::vector<std::vector<Var>> &equivVar) {

  constexpr int max_depth = 3;
  constexpr double max_time = 3;
  int tested = 0;
  struct hashFunction {
    size_t operator()(const std::vector<Var> &myVector) const {
      std::hash<bool> hasher;
      size_t answer = 0;
      for (int i : myVector) {
        answer ^= hasher(i) + 0x9e3779b9 + (answer << 6) + (answer >> 2);
      }
      return answer;
    }
  };

  struct EquivClass {
    std::vector<Var> rep;
    std::vector<Var> vars;
    int nb_proj;
    int nb_nproj;
    bool bad;
    float gain() const { return float(nb_nproj + nb_proj) / rep.size(); }
  };
  auto count_proj = [&](std::vector<Var> &vars) {
    int i = 0;
    for (auto v : vars) {
      i += m_specs->isSelected(v);
    }
    return i;
  };
  auto create_equiv = [&](std::vector<Var> &&vars, std::vector<Var> &&rep,
                          bool bad = false) {
    int nb_proj = count_proj(vars);
    int nb_nproj = vars.size() - nb_proj;
    return EquivClass{std::move(rep), std::move(vars), nb_proj, nb_nproj, bad};
  };
  auto cmp_equiv_class = [&](const EquivClass &lhs, const EquivClass &rhs) {
    if (lhs.nb_nproj == rhs.nb_nproj) {
      return lhs.nb_proj < rhs.nb_proj;
    }
    return lhs.nb_nproj > rhs.nb_nproj;
  };

  auto merge_equiv = [&]() {
    std::fill(m_markedVar.begin(), m_markedVar.end(), false);
    std::fill(m_flagVar.begin(), m_flagVar.end(), false);
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
  };
  auto merge_vectors = [](std::vector<int> &a, std::vector<int> &b) {
    std::vector<int> c = a;
    c.insert(c.end(), b.begin(), b.end());
    std::sort(c.begin(), c.end());
    c.erase(std::unique(c.begin(), c.end()), c.end());
    return c;
  };

  std::vector<Var> proj_vars(vars.begin(), vars.begin() + count_proj(vars));

  std::sort(proj_vars.begin(), proj_vars.end(), [&](Var a, Var b) {
    return m_specs->getNbOccurrence(Lit::makeLitTrue(a)) +
               m_specs->getNbOccurrence(Lit::makeLitFalse(a)) >
           m_specs->getNbOccurrence(Lit::makeLitTrue(b)) +
               m_specs->getNbOccurrence(Lit::makeLitFalse(b));
  });
  std::vector<EquivClass> current_eq;
  std::vector<std::vector<bool>> graph_index(
      m_specs->getNbVariable() + 1,
      std::vector<bool>(m_specs->getNbVariable() + 1, false));
  std::vector<std::vector<Var>> graph(m_specs->getNbVariable() + 1);
  std::vector<unsigned int> clauses;
  m_specs->getCurrentClauses(clauses, vars);
  for (auto cidx : clauses) {
    auto cl = m_specs->getClause(cidx);
    if (cl.size() > 2)
      continue;
    for (int i = 0; i < cl.size(); i++) {
      if (!m_specs->isSelected(cl[i].var()))
        continue;
      for (int j = i + 1; j < cl.size(); j++) {
        if (!m_specs->isSelected(cl[j].var()))
          continue;
        if (!graph[cl[i].var()][cl[j].var()]) {
          graph[cl[i].var()].push_back(cl[j].var());
          graph[cl[i].var()][cl[j].var()] = true;
        }
      }
    }
  }

  std::vector<std::vector<int>> hits(
      max_depth + 1, std::vector<int>(m_specs->getNbVariable() + 1));
  std::vector<int> sat_cnt(m_specs->getNbVariable() + 1);
  auto clear_cnts = [&]() {
    for (auto &k : hits) {
      for (auto &h : k) {
        h = 0;
      }
    }
    for (auto &k : sat_cnt) {
      k = 0;
    }
  };
  for (Var v : proj_vars) {
    auto &bin_neigh = graph[v];
    if (bin_neigh.size() < max_depth - 1) {
      std::vector<Var> sel({v});
      for (auto l : bin_neigh) {
        sel.push_back(l);
      }
      clear_cnts();
      search_units(s,vars,sel,0,hits,sat_cnt);
      for(int i = 1;i<= max_depth;i++){


      }
    } else {
    }
  }
  {
    // Create all D=1 classes
    for (auto &v : vars)
      m_flagVar[v] = true;
    for (auto &v : proj_vars) {
      assert((unsigned)v < m_markedVar.size());
      if (m_markedVar[v] || s.varIsAssigned(v))
        continue;
      std::vector<Var> eqv;
      if (interCollectUnit(s, v, eqv, m_flagVar)) {
        assert(eqv.size() > 0);
        if (eqv.size() == 1)
          continue;
        for (auto &vv : eqv) {
          m_markedVar[vv] = true;
        }
        current_eq.push_back(create_equiv(std::move(eqv), {v}));
      }
    }
  }
  if (proj_vars.size() > 16) {
    auto start = std::chrono::steady_clock::now();
    std::unordered_set<std::vector<Var>, hashFunction> already_seen;
    for (int d = 0; d < max_depth - 1; d++) {
      for (auto v : proj_vars) {
        auto finish = std::chrono::steady_clock::now();
        double elapsed_seconds =
            std::chrono::duration_cast<std::chrono::duration<double>>(finish -
                                                                      start)
                .count();
        if (elapsed_seconds > max_time) {
          goto END;
        }
      }
    }
  }
END:
  std::sort(current_eq.begin(), current_eq.end(),
            [](const EquivClass &lhs, const EquivClass &rhs) {
              return lhs.gain() < rhs.gain();
            });

  std::fill(m_markedVar.begin(), m_markedVar.end(), false);
  int nb_proj = 0;
  for (int i = current_eq.size() - 1; i >= 0; i--) {
    bool remove = current_eq[i].vars.size() == 1;
    for (auto r : current_eq[i].vars) {
      if (m_markedVar[r]) {
        remove = true;
      }
    }
    if (remove) {
      std::swap(current_eq[i], current_eq.back());
      current_eq.pop_back();
    } else {
      for (Var v : current_eq[i].vars) {
        m_markedVar[v] = true;
        nb_proj++;
      }
    }
  }
  std::fill(m_markedVar.begin(), m_markedVar.end(), false);
  std::cout << "Classes\n";
  for (auto &e : current_eq) {
    std::sort(e.vars.begin(), e.vars.end());
    std::cout << e.rep << "=>" << e.vars << std::endl;
    equivVar.push_back(std::move(e.vars));
  }

  // merge_equiv();
  std::cout << "Merged:" << std::endl;
  for (auto &r : equivVar) {
    std::cout << r << std::endl;
  }

  std::cout << "Tested " << tested << std::endl;
  return;

  /*
    auto get_cut = [&]() {
      std::vector<Var> equivClass(m_specs->getNbVariable() + 1);
      std::vector<unsigned> var2class(m_specs->getNbVariable() + 1, -1);
      for (auto &v : vars) {
        assert(equivClass.size() >= (unsigned)v);
        equivClass[v] = v;
      }
      for (auto &c : equivVar) {
        Var vi = c.back();
        for (auto &v : c)
          equivClass[v] = vi;
      }
      std::vector<Var> considered;
      HyperGraph h;
      m_hex->constructHyperGraph(*m_specs, vars, equivClass, equivVar, true,
                                 considered, h);
      std::vector<int> partition;
      m_part->computePartition(h, PartitionerManager::QUALITY, partition);
      std::vector<unsigned> cut;
      for (auto &e : h) {
        int p = partition[e.getId()];
        bool in_cut = false;
        for (auto v : e) {
          if (p != partition[v]) {
            in_cut = true;
          }
        }
        if (in_cut) {
          cut.push_back(e.getId());
        }
      }
      for (int i = cut.size() - 1; i >= 0; i--) {
        if (var2class[cut[i]] == -1) {
          continue;
        }
        for (auto e : equivVar[var2class[cut[i]]]) {
          cut.push_back(e);
        }
      }
      for (int i = cut.size() - 1; i >= 0; i--) {
        if (m_specs->isSelected(cut[i])) {
          std::swap(cut[i], cut.back());
          cut.pop_back();
        }
      }
      return cut;
    };
  */

} // searchEquiv

} // namespace d4
