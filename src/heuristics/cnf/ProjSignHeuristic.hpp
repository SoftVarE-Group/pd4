#pragma once

#include <ostream>

#include "PhaseSelectorManager.hpp"
namespace d4 {
class ProjSignHeuristic {
public:
  Var select(std::vector<d4::Var> &vars, std::vector<bool> &priority);

  ProjSignHeuristic() : m_om(nullptr) {}
  ProjSignHeuristic(SpecManagerCnf *om) : m_om(om) {}

private:
  std::vector<unsigned int> m_clauses;
  d4::SpecManagerCnf *m_om;
};
} // namespace d4
