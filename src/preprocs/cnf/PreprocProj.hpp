#include "util/TestSolver.hpp"
#include "util/lib_sharpsat_td/bitset.hpp"
#include "util/lib_sharpsat_td/subsumer.hpp"
#include "util/lib_sharpsat_td/utils.hpp"
#include <algorithm>

#pragma once

#include <boost/program_options.hpp>
#include <vector>

#include "../PreprocManager.hpp"
#include "3rdParty/glucose-3.0/core/SolverTypes.h"
#include "src/problem/ProblemTypes.hpp"
#include "src/solvers/WrapperSolver.hpp"

namespace d4 {
namespace po = boost::program_options;

class PreprocProj : public PreprocManager {
private:
  struct Config {
      int reps = 20;
      double max_time = 1.0;
  } config;
  struct Stats{
      int elim_mono = 0;
      int elim_stren =0;
      int elim_adj =0;
      int elim_adj_var =0;

  } stats;
  std::vector<std::vector<Glucose::Lit>> clauses, learnts;
  int vars, pvars, freevars;
  std::vector<Glucose::Lit> assignedLits;
  std::vector<Glucose::lbool> assigns;
  std::vector<bool> ispvars;
  std::vector<Glucose::Lit>
      gmap; // variable renaming (original (Var) -> simplified (Lit))
  std::vector<Glucose::Lit> fixedLits; // fixed Literals
  std::vector<Glucose::Var>
      definedVars; // variables defined by the other variables
  std::vector<std::vector<Glucose::Lit>>
      freeLitClasses; // equivalence classes of free literals
  bool unsat;
  WrapperSolver* ws;

public:
  PreprocProj(po::variables_map &vm, std::ostream &out);
  ~PreprocProj();
  bool sat_fle();
  void strenthen();
  bool elim_mono();
  void compact(const Glucose::vec<Glucose::lbool> &assigns,
               const std::vector<Glucose::Var> &elimvars = {});
  void subsume();
  void compact_clauses(const Glucose::vec<Glucose::lbool> &assigns,
                       std::vector<std::vector<Glucose::Lit>> &cls,
                       std::vector<bool> &occurred, int &varnum);
  void rewrite_claues(std::vector<std::vector<Glucose::Lit>> &cls,
                      const std::vector<Glucose::Lit> &map);
  void rewrite_claues(std::vector<std::vector<Glucose::Lit>> &cls,
                      const std::vector<Glucose::Var> &map);
  bool add_clause(std::vector<Glucose::Lit> &clause);
  void merge_equiv();
  Glucose::lbool value(Glucose::Lit);
  void simplify();
  virtual ProblemManager *run(ProblemManager *pin,
                              LastBreathPreproc &lastBreath) override;
};
} // namespace d4
