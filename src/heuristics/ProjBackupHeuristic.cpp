#include "ProjBackupHeuristic.hpp"
#include "src/heuristics/cnf/ProjBackupHeuristicHypergraph.hpp"
#include <string>

#include "src/heuristics/cnf/ProjBackupHeuristicComponents.hpp"

namespace d4 {

std::unique_ptr<ProjBackupHeuristic>
ProjBackupHeuristic::make(po::variables_map &vm, SpecManager &p,
                          WrapperSolver &s, std::ostream &out) {
  std::string meth = vm["proj-backup"].as<std::string>();
  if (meth == "hypergraph") {
    return std::make_unique<ProjBackupHeuristicHypergraph>(vm, p, s, out);

  } else if (meth == "component") {
    return std::make_unique<ProjBackupHeuristicComponents>(vm, p, s, out);

  } else if (meth == "none") {
    return std::make_unique<ProjBackupHeuristic>();
  }
  throw std::runtime_error("unknown proj-backup method");
}
bool ProjBackupHeuristic::computeCutSetDyn(ProjVars &component,
                                           std::vector<Var> &cutset) {
  return false;
}
} // namespace d4
