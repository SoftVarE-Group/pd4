
#include <ostream>

#include "ProjSignHeuristic.hpp"
#include <memory>
namespace d4 {
int sign_matches[1 << 14][1 << 10][2];
int priority_map[1 << 14];
Var ProjSignHeuristic::select(std::vector<d4::Var> &vars,
                              std::vector<bool> &priority) {

  std::memset(sign_matches, 0, sizeof(sign_matches));
  int priority_cnt = 0;
  for (Var v : vars) {
    if (priority[v]) {
      priority_map[v] = priority_cnt;
      priority_cnt++;
    }
  }
  if (priority_cnt > 1 << 10) {
    return var_Undef;
  }
  m_om->getCurrentClauses(m_clauses, vars);

  std::vector<Lit> priority_in_clause;

  for (auto ci : m_clauses) {
    priority_in_clause.clear();
    for (auto l : m_om->getClause(ci)) {
      if (m_om->litIsAssigned(l))
        continue;
      if (priority[l.var()]) {
        priority_in_clause.push_back(l);
      }
    }
    if (priority_in_clause.empty()) {
      continue;
    }

    for (auto l : m_om->getClause(ci)) {
      if (m_om->litIsAssigned(l) || !m_om->isSelected(l.var()))
        continue;
      for (auto k : priority_in_clause) {
        int s = k.sign() ^ l.sign();
        sign_matches[l.var()][priority_map[k.var()]][s]++;
      }
    }
  }
  int direct_matches = 0;
  int direct_matches_var = var_Undef;
  double indirect_macthes_ration = 1.0;
  int indirect_macthes_var = var_Undef;

  for (auto v : vars) {
    if (m_om->varIsAssigned(v) || !m_om->isSelected(v))
      continue;
    for (int i = 0; i < priority_cnt; i++) {
      int a = sign_matches[v][i][0];
      int b = sign_matches[v][i][1];
      if (a == 0 && b == 0) {
        continue;
      } else if (a != 0 && b == 0) {
        if (a > direct_matches) {
          direct_matches = a;
          direct_matches_var = v;
        }
      } else if (a == 0 && b != 0) {
        if (b > direct_matches) {
          direct_matches = b;
          direct_matches_var = v;
        }
      } else {
        double ratio = abs(a - b) / double(a + b);
        if (ratio > indirect_macthes_ration) {
          indirect_macthes_ration = ratio;
          indirect_macthes_var = v;
        }
      }
    }
  }
  if (direct_matches_var == var_Undef) {
    return var_Undef;

  } else {
    std::cout << "Found direct match" << std::endl;
    return direct_matches_var;
  }
}

} // namespace d4
